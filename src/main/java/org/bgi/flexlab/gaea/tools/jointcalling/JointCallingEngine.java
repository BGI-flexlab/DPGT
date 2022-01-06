package org.bgi.flexlab.gaea.tools.jointcalling;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.annotator.RMSAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.annotator.RankSumTest;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.*;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCallingOptions;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.seqdoop.hadoop_bam.LazyParsingGenotypesContext;
import org.seqdoop.hadoop_bam.LazyVCFGenotypesContext.HeaderDataCache;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class JointCallingEngine {
	
	private static final String GVCF_BLOCK = "GVCFBlock";
	private final List<String> infoFieldAnnotationKeyNamesToRemove = new ArrayList<>();
	
	private final boolean INCLUDE_NON_VARIANTS;

	private final boolean uniquifySamples;

	private final TreeMap<Integer, ArrayList<VariantContext>> variantsForSample;

	private VariantContext currentContext = null;
	private static Set<VCFHeaderLine> gvcfHeaderMetaInfo;
	private int max_position = -1;

	// the genotyping engine
	private final UnifiedGenotypingEngine genotypingEngine;
	// the annotation engine
	private final VariantAnnotatorEngine annotationEngine;

	private final GenomeLocationParser parser;

	protected List<String> annotationsToUse = new ArrayList<>();

	protected List<String> annotationGroupsToUse = new ArrayList<>(
			Collections.singletonList(StandardAnnotation.class.getSimpleName()));
	private final VCFHeader vcfHeader;
	private final HashMap<Integer,Set<String>> nameHeaders= new HashMap<>();
	final Set<String> infoHeaderAltAllelesLineNames = new LinkedHashSet<>();
	private final Set<Integer> mapMergedSamples= new TreeSet<>();
	private int sampleSize = 0;
	public JointCallingEngine(JointCallingOptions options, GenomeLocationParser parser,
			MultipleVCFHeaderForJointCalling multiHeaders) {
		variantsForSample = new TreeMap<>();
		this.INCLUDE_NON_VARIANTS = options.INCLUDE_NON_VARIANT;
		this.uniquifySamples = options.isUniquifySamples();
		this.parser = parser;

		annotationEngine = new VariantAnnotatorEngine(annotationGroupsToUse, annotationsToUse,
				Collections.emptyList());
		annotationEngine.initializeDBs(options.getDBSnp() != null);
		
		for ( final InfoFieldAnnotation annotation :  annotationEngine.getRequestedInfoAnnotations() ) {
            if ( annotation instanceof RankSumTest || annotation instanceof RMSAnnotation ) {
                final List<String> keyNames = annotation.getKeyNames();
                if ( !keyNames.isEmpty() ) {
                    infoFieldAnnotationKeyNamesToRemove.add(keyNames.get(0));
                }
            }
        }
		Set<String> sampleNames = getSampleList(multiHeaders.getMergeHeader());

		genotypingEngine = new UnifiedGenotypingEngine(sampleNames.size(), options, this.parser);

		// take care of the VCF headers
		final Set<VCFHeaderLine> headerLines = new HashSet<>();

		headerLines.addAll(multiHeaders.getMergeHeader().getMetaDataInInputOrder());

		headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions());
		headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

		// add headers for annotations added by this tool
		headerLines.add(GaeaVcfHeaderLines.getInfoLine(GaeaVCFConstants.MLE_ALLELE_COUNT_KEY));
		headerLines.add(GaeaVcfHeaderLines.getInfoLine(GaeaVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
		headerLines.add(GaeaVcfHeaderLines.getFormatLine(GaeaVCFConstants.REFERENCE_GENOTYPE_QUALITY));
		headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));
		
		if ( INCLUDE_NON_VARIANTS ) {
            // Save INFO header names that require alt alleles
            for ( final VCFHeaderLine headerLine : headerLines ) {
                if (headerLine instanceof VCFInfoHeaderLine ) {
                    if (((VCFInfoHeaderLine) headerLine).getCountType() == VCFHeaderLineCount.A) {
                        infoHeaderAltAllelesLineNames.add(((VCFInfoHeaderLine) headerLine).getID());
                    }
                }
            }
        }
		
		if (options.getDBSnp() != null)
			VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);

		vcfHeader = new VCFHeader(headerLines, sampleNames);

		this.sampleSize = multiHeaders.getHeaderSize();

		Set<String> sampleNamesHashSet = new HashSet<>(multiHeaders.getMergeHeader().getSampleNamesInOrder());
		annotationEngine.invokeAnnotationInitializationMethods(headerLines, sampleNamesHashSet);
		gvcfHeaderMetaInfo=multiHeaders.getMergeHeader().getMetaDataInInputOrder();
		GvcfMathUtils.resetRandomGenerator();

	}

	public JointCallingEngine(JointCallingOptions options, GenomeLocationParser parser, VCFHeader vcfheader,
							  MultipleVCFHeaderForJointCalling multiHeaders, ArrayList<ArrayList<String>> multiMapSamples, Configuration conf) throws IllegalArgumentException, IOException{
		variantsForSample = new TreeMap<>();
		this.INCLUDE_NON_VARIANTS = options.INCLUDE_NON_VARIANT;
		this.uniquifySamples = options.isUniquifySamples();
		this.parser = parser;

		annotationEngine = new VariantAnnotatorEngine(annotationGroupsToUse, annotationsToUse,
				Collections.emptyList());
		annotationEngine.initializeDBs(options.getDBSnp() != null);
		
		for ( final InfoFieldAnnotation annotation :  annotationEngine.getRequestedInfoAnnotations() ) {
            if ( annotation instanceof RankSumTest || annotation instanceof RMSAnnotation ) {
                final List<String> keyNames = annotation.getKeyNames();
                if ( !keyNames.isEmpty() ) {
                    infoFieldAnnotationKeyNamesToRemove.add(keyNames.get(0));
                }
            }
        }
		Set<String> sampleNames = getSampleList(vcfheader);

		genotypingEngine = new UnifiedGenotypingEngine(sampleNames.size(), options, this.parser);

		// take care of the VCF headers

		final Set<VCFHeaderLine> headerLines = new HashSet<>(multiHeaders.getMergeHeader().getMetaDataInInputOrder());

		headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().contains(GVCF_BLOCK));
		for(VCFHeaderLine hl:annotationEngine.getVCFAnnotationDescriptions()){
            headerLines.removeIf(vcfHeaderLine -> idEquals(vcfHeaderLine,hl));
            headerLines.add(hl);
        }
        for(VCFHeaderLine hl:genotypingEngine.getAppropriateVCFInfoHeaders()){
            headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().equals(hl.getKey()));
            headerLines.add(hl);
        }
		headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().contains(GVCF_BLOCK));
		headerLines.removeIf(vcfHeaderLine -> (vcfHeaderLine.getKey().equals("INFO") && getInfoID(vcfHeaderLine).equals(GaeaVCFConstants.MLE_ALLELE_COUNT_KEY)));
		headerLines.removeIf(vcfHeaderLine -> (vcfHeaderLine.getKey().equals("INFO") && getInfoID(vcfHeaderLine).equals(GaeaVCFConstants.MLE_ALLELE_FREQUENCY_KEY)));
		headerLines.removeIf(vcfHeaderLine -> (vcfHeaderLine.getKey().equals("INFO") && getInfoID(vcfHeaderLine).equals(GaeaVCFConstants.REFERENCE_GENOTYPE_QUALITY)));
		headerLines.removeIf(vcfHeaderLine -> (vcfHeaderLine.getKey().equals("INFO") && getInfoID(vcfHeaderLine).equals(VCFConstants.DEPTH_KEY)));
		headerLines.add(GaeaVcfHeaderLines.getInfoLine(GaeaVCFConstants.MLE_ALLELE_COUNT_KEY));
		headerLines.add(GaeaVcfHeaderLines.getInfoLine(GaeaVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
		headerLines.add(GaeaVcfHeaderLines.getFormatLine(GaeaVCFConstants.REFERENCE_GENOTYPE_QUALITY));
		headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));
		if ( INCLUDE_NON_VARIANTS ) {
            // Save INFO header names that require alt alleles
            for ( final VCFHeaderLine headerLine : headerLines ) {
                if (headerLine instanceof VCFInfoHeaderLine ) {
                    if (((VCFInfoHeaderLine) headerLine).getCountType() == VCFHeaderLineCount.A) {
                        infoHeaderAltAllelesLineNames.add(((VCFInfoHeaderLine) headerLine).getID());
                    }
                }
            }
        }
		
		if (options.getDBSnp() != null)
			VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
		vcfHeader = new VCFHeader(headerLines, sampleNames);
//		
		ArrayList<String> totalSamples= new ArrayList<>();
		HashMap<String, Integer> sampleIndex=new HashMap<>();
		String sampleIndexFile=options.getOutDir()+"/vcfheaderinfo";
		Path sampleIndexFilePath=new Path(sampleIndexFile);
		FileSystem fs=sampleIndexFilePath.getFileSystem(conf);
		BufferedReader reader;
		if(sampleIndexFile.startsWith("hdfs:")){
			reader = new BufferedReader(new InputStreamReader(fs.open(sampleIndexFilePath)));
		}else {
			if(sampleIndexFile.startsWith("file://")){
				sampleIndexFile=sampleIndexFile.substring(7);
			}
			reader = new BufferedReader(new FileReader(sampleIndexFile));
		}
		String line;
		Logger logger=LoggerFactory.getLogger(JointCallingEngine.class);
		int totalSampleSize=0;
		while((line=reader.readLine())!=null) {
			totalSampleSize++;
			String[] eles=line.split("\t");
			if(eles.length!=3) {
				logger.error("vcfheaderinfo file format error");
			}
			String name;
			if(eles[2].endsWith(",")) {
				name=eles[2].substring(0,eles[2].length()-1);
			}else {
				name=eles[2];
			}
			sampleIndex.put(name, Integer.parseInt(eles[1]));
			Set<String> samples = new LinkedHashSet<>();
			String[] sample = eles[2].split(",");
			Collections.addAll(samples, sample);
			nameHeaders.put(Integer.parseInt(eles[1]), samples);
		}
		reader.close();
		//create nameHeaders

		File rawInput=new File(options.getInputList());
		String sortedInputGvcfList=rawInput.getParent()+"/sorted."+rawInput.getName();
		if(sortedInputGvcfList.startsWith("file:")){
			sortedInputGvcfList=sortedInputGvcfList.substring(5);
		}
		BufferedReader inputListReader;
		File sparkOrMapreducePosition=new File(sortedInputGvcfList);
		if(sparkOrMapreducePosition.exists()) {
			if (!sortedInputGvcfList.startsWith("file://")) {
				sortedInputGvcfList = "file://" + sortedInputGvcfList;
			}
			Path inputListPath = new Path(sortedInputGvcfList);
			inputListReader = new BufferedReader(new InputStreamReader(inputListPath.getFileSystem(conf).open(inputListPath)));
		}else{
			sortedInputGvcfList=options.getOutDir()+"/sorted."+rawInput.getName();
			inputListReader=new BufferedReader(new FileReader(sortedInputGvcfList));
		}

		inputListReader.close();
		for(ArrayList<String> mapSamples:multiMapSamples) {
			Set<String> samples = new LinkedHashSet<>(mapSamples);
			Integer mergedSampleIndex=0;
			int mapSamplesIndex=0;
            for(String sm:mapSamples) {
            	totalSamples.add(sm);
            	if(mapSamplesIndex==0) {
					mergedSampleIndex = sampleIndex.get(sm);
				}
				mapSamplesIndex++;
            }
			mergedSampleIndex+=totalSampleSize;
            nameHeaders.put(mergedSampleIndex, samples);
            mapMergedSamples.add(mergedSampleIndex);
		}
		Set<String> sampleNamesHashSet = new HashSet<>(totalSamples);
		annotationEngine.invokeAnnotationInitializationMethods(headerLines, sampleNamesHashSet);
		gvcfHeaderMetaInfo=multiHeaders.getMergeHeader().getMetaDataInInputOrder();
		GvcfMathUtils.resetRandomGenerator();
	}

	private String getInfoID(VCFHeaderLine vcfHeaderLine) {
		String pattern=".*<ID=(.*?),.*";
		Pattern r=Pattern.compile(pattern);
		Matcher h1m=r.matcher(vcfHeaderLine.toString());
		if(h1m.find()){
			return h1m.group(1);
		}else{
			return "";
		}
	}

	private boolean idEquals(VCFHeaderLine vcfHeaderLine, VCFHeaderLine hl) {
		if(vcfHeaderLine.getKey().equals(hl.getKey())){
			String h1ID=getInfoID(vcfHeaderLine);
			String h2ID=getInfoID(hl);
			if(!Objects.equals(h1ID, "") && !Objects.equals(h2ID, "")){
				return h1ID.equals(h2ID);
			}else{
				return false;
			}
		}else{
			return false;
		}
    }

    public static Set<String> getSampleList(Set<VCFHeader> headers) {
		Set<String> samples = new TreeSet<>();
		for (VCFHeader header : headers) {
			for (String sample : header.getGenotypeSamples()) {
				samples.add(GaeaGvcfVariantContextUtils.mergedSampleName(null, sample, false));
			}
		}

		return samples;
	}

	public Set<String> getSampleList(VCFHeader header) {
		Set<String> samples = new TreeSet<>();
		for (String sample : header.getGenotypeSamples()) {
			samples.add(GaeaGvcfVariantContextUtils.mergedSampleName(null, sample, false));
		}

		return samples;
	}

	public void init(ArrayList<VariantContext> dbsnps) {
		annotationEngine.initializeDBs(dbsnps, parser);
	}

	public void purgeOutOfScopeRecords(GenomeLocation location) {
		for (Integer sample : variantsForSample.keySet()) {
			variantsForSample.get(sample).removeIf(context -> context.getEnd() < location.getStart());
		}
	}

	public void lazyLoad(Iterator<VariantContextWritable> iterator, GenomeLocation location) {
		int curr = location.getStart();

		if (curr <= max_position)
			purgeOutOfScopeRecords(location);
		else {
			for (Integer sample : variantsForSample.keySet()) {
				variantsForSample.get(sample).clear();
			}
			max_position = -1;
		}

		if (currentContext == null) {
			if (iterator.hasNext()) {
				currentContext = iterator.next().get();
			}
		}
		while (currentContext != null) {
			if (currentContext.getStart() > curr)
				break;
			if (currentContext.getStart() <= curr && currentContext.getEnd() >= curr) {
				GenotypesContext gc = currentContext.getGenotypes();
				String sampleName = currentContext.getAttributeAsString("SM", null);
				if (sampleName == null)
					throw new UserException("Not contains SM attribute");
				
				int sampleID = Integer.parseInt(sampleName);
				
				if (gc instanceof LazyParsingGenotypesContext) {
					VCFHeader header = new VCFHeader(gvcfHeaderMetaInfo, nameHeaders.get(sampleID));
					HeaderDataCache datacache = new HeaderDataCache();
					datacache.setHeader(header);
					((LazyParsingGenotypesContext) gc).getParser().setHeaderDataCache(datacache);
				}
				
				if (variantsForSample.containsKey(sampleID)) {
					variantsForSample.get(sampleID).add(currentContext);
				} else {
					ArrayList<VariantContext> list = new ArrayList<>();
					list.add(currentContext);
					variantsForSample.put(sampleID, list);
				}

				if (max_position < currentContext.getEnd())
					max_position = currentContext.getEnd();
			}

			if (iterator.hasNext()) {
				currentContext = iterator.next().get();
			} else {
				currentContext = null;
				max_position = -1;
			}
		}
	}
	public void lazyLoadForSpark(Iterator<VariantContext> iterator, GenomeLocation location) {
		int curr = location.getStart();

		if (curr <= max_position)
			purgeOutOfScopeRecords(location);
		else {
			for (Integer sample : variantsForSample.keySet()) {
				variantsForSample.get(sample).clear();
			}
			max_position = -1;
		}

		if (currentContext == null) {
			if (iterator.hasNext()) {
				currentContext = iterator.next();
			}
		}
		while (currentContext != null) {
			if (currentContext.getStart() > curr)
				break;
			if (currentContext.getStart() <= curr && currentContext.getEnd() >= curr) {
				GenotypesContext gc = currentContext.getGenotypes();
				String sampleName = currentContext.getAttributeAsString("SM", null);
				if (sampleName == null) {
					throw new UserException("Not contains SM attribute");
				}
				int sampleID = Integer.parseInt(sampleName);

				if (gc instanceof LazyParsingGenotypesContext) {
					VCFHeader header = new VCFHeader(gvcfHeaderMetaInfo, nameHeaders.get(sampleID));
					HeaderDataCache datacache = new HeaderDataCache();
					datacache.setHeader(header);
					((LazyParsingGenotypesContext) gc).getParser().setHeaderDataCache(datacache);
				}

				if (variantsForSample.containsKey(sampleID)) {
					variantsForSample.get(sampleID).add(currentContext);
				} else {
					ArrayList<VariantContext> list = new ArrayList<>();
					list.add(currentContext);
					variantsForSample.put(sampleID, list);
				}

				if (max_position < currentContext.getEnd())
					max_position = currentContext.getEnd();
			}

			if (iterator.hasNext()) {
				currentContext = iterator.next();
			} else {
				currentContext = null;
				max_position = -1;
			}
		}
	}

	private VariantContext getValues(int sample,GenomeLocation loc,boolean requireStartHere){
		if(variantsForSample.containsKey(sample) && variantsForSample.get(sample).size() > 0){
			for(VariantContext vc : variantsForSample.get(sample)){
				if ( ! requireStartHere || vc.getStart() == loc.getStart()){
					return vc;
				}
			}
		}
		
		return null;
	}
	
	private VariantContext getValues(int sample,GenomeLocation loc){
		VariantContext vc = getValues(sample,loc,true);
		if(vc == null)
			vc = getValues(sample,loc,false);
		
		return vc;
	}

	private List<VariantContext> getValues(GenomeLocation loc) {
		List<VariantContext> list = new ArrayList<>();
		
		int i;
		for(i = 0 ; i < sampleSize ; i++) {
			VariantContext vc = getValues(i,loc);
			if(vc != null)
				list.add(vc);
		}
		for(Integer index:mapMergedSamples) {
			VariantContext vc = getValues(index,loc);
			if(vc != null)
				list.add(vc);
		}
		return list;
	}

	public VariantContext variantCallingReduce(Iterator<VariantContextWritable> iterator, GenomeLocation location,
										 ChromosomeInformationShare ref) {
		if (location.getStart() != location.getStop())
			throw new UserException("location must length is 1!");
		lazyLoad(iterator, location);

		final List<VariantContext> vcsAtThisLocus = getValues(location);
		final Byte refBase = INCLUDE_NON_VARIANTS ? (byte) ref.getBase(location.getStart() - 1) : null;
		final boolean removeNonRefSymbolicAllele = !INCLUDE_NON_VARIANTS;
		final VariantContext combinedVC = ReferenceConfidenceVariantContextMerger.reduceMerge(vcsAtThisLocus, location,
				refBase, removeNonRefSymbolicAllele, uniquifySamples, annotationEngine);
		return combinedVC == null ? null : regenotypeVC(new RefMetaDataTracker(location), ref, combinedVC);
	}
	public VariantContext variantCallingForSpark(Iterator<VariantContext> iterator, GenomeLocation location,
												 ChromosomeInformationShare ref) {
		if (location.getStart() != location.getStop())
			throw new UserException("location must length is 1!");
		lazyLoadForSpark(iterator, location);

		final List<VariantContext> vcsAtThisLocus = getValues(location);
		final Byte refBase = INCLUDE_NON_VARIANTS ? (byte) ref.getBase(location.getStart() - 1) : null;
		final boolean removeNonRefSymbolicAllele = !INCLUDE_NON_VARIANTS;
		final VariantContext combinedVC = ReferenceConfidenceVariantContextMerger.reduceMerge(vcsAtThisLocus, location,
				refBase, removeNonRefSymbolicAllele, uniquifySamples, annotationEngine);
		return combinedVC == null ? null : regenotypeVC(new RefMetaDataTracker(location), ref, combinedVC);
	}

	protected VariantContext regenotypeVC(final RefMetaDataTracker tracker, final ChromosomeInformationShare ref,
			final VariantContext originalVC) {
		if (originalVC == null) {
			throw new IllegalArgumentException("originalVC cannot be null");
		} else if (!isProperlyPolymorphic(originalVC) && !INCLUDE_NON_VARIANTS) {
			return null;
		}

		VariantContext result = originalVC;
		// don't need to calculate quals for sites with no data whatsoever
		if (result.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > 0) {
			result = genotypingEngine.calculateGenotypes(originalVC);
		}
		if (result == null || (!isProperlyPolymorphic(result) && !INCLUDE_NON_VARIANTS)) {
			return null;
		}
		result = addGenotypingAnnotations(originalVC.getAttributes(), result);
		// At this point we should already have DP and AD annotated
		result = annotationEngine.finalizeAnnotations(result, originalVC);
		// do trimming after allele-specific annotation reduction or the mapping
		// is difficult
		result = GaeaGvcfVariantContextUtils.reverseTrimAlleles(result);
		// Re-annotate and fix/remove some of the original annotations.
		// Note that the order of these actions matters and is different for
		// polymorphic and monomorphic sites.
		// For polymorphic sites we need to make sure e.g. the SB tag is sent to
		// the annotation engine and then removed later.
		// For monomorphic sites we need to make sure e.g. the hom ref genotypes
		// are created and only then are passed to the annotation engine.
		// We could theoretically make 2 passes to re-create the genotypes, but
		// that gets extremely expensive with large sample sizes.
		if (result.isPolymorphicInSamples()) {
			result = annotationEngine.annotateContext(tracker, ref, result);
			result = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, false)).make();
		} else if (INCLUDE_NON_VARIANTS) {
			result = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, true)).make();
			result = annotationEngine.annotateContext(tracker, ref, result);
			result = removeNonRefAlleles(result);
		} else {
			return null;
		}
		result = removeInfoAnnotationsIfNoAltAllele(result);
		return result;
	}
	
	 /**
     * Remove INFO field annotations if no alternate alleles
     *
     * @param vc    the variant context
     * @return variant context with the INFO field annotations removed if no alternate alleles
    */
    private VariantContext removeInfoAnnotationsIfNoAltAllele(final VariantContext vc)  {

        // If no alt alleles, remove any RankSumTest or RMSAnnotation attribute
        if ( vc.getAlternateAlleles().isEmpty() ) {
            final VariantContextBuilder builder = new VariantContextBuilder(vc);

            for ( final String annotation : infoFieldAnnotationKeyNamesToRemove ) {
                builder.rmAttribute(annotation);
            }
            return builder.make();
        } else {
            return vc;
        }
    }
	
	/**
     * Remove NON-REF alleles from the variant context
     *
     * @param vc   the variant context
     * @return variant context with the NON-REF alleles removed if multiallelic or replaced with NO-CALL alleles if biallelic
     */
    private VariantContext removeNonRefAlleles(final VariantContext vc) {

        // If NON_REF is the only alt allele, ignore this site
        final List<Allele> newAlleles = new ArrayList<>();
        // Only keep alleles that are not NON-REF
        for ( final Allele allele : vc.getAlleles() ) {
            if ( !allele.equals(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE) ) {
                newAlleles.add(allele);
            }
        }

        // If no alt allele, then remove INFO fields that require alt alleles
        if ( newAlleles.size() == 1 ) {
            final VariantContextBuilder builder = new VariantContextBuilder(vc).alleles(newAlleles);
            for ( final String name : infoHeaderAltAllelesLineNames ) {
                builder.rmAttributes(Collections.singletonList(name));
            }
            return builder.make();
        } else {
            return vc;
        }
    }

	private boolean isProperlyPolymorphic(final VariantContext vc) {
		// obvious cases
		if (vc == null || vc.getAlternateAlleles().isEmpty()) {
			return false;
		} else if (vc.isBiallelic()) {
			return !(vc.getAlternateAllele(0).equals(Allele.SPAN_DEL)
					|| vc.getAlternateAllele(0).equals(GaeaVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED)
					|| vc.isSymbolic());
		} else {
			return true;
		}
	}

	private VariantContext addGenotypingAnnotations(final Map<String, Object> originalAttributes,
			final VariantContext newVC) {
		// we want to carry forward the attributes from the original VC but make
		// sure to add the MLE-based annotations
		final Map<String, Object> attrs = new HashMap<>(originalAttributes);
		attrs.put(GaeaVCFConstants.MLE_ALLELE_COUNT_KEY, newVC.getAttribute(GaeaVCFConstants.MLE_ALLELE_COUNT_KEY));
		attrs.put(GaeaVCFConstants.MLE_ALLELE_FREQUENCY_KEY,
				newVC.getAttribute(GaeaVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
		if (newVC.hasAttribute(GaeaVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY))
			attrs.put(GaeaVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY,
					newVC.getAttribute(GaeaVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));
		if (newVC.hasAttribute(GaeaVCFConstants.AS_QUAL_KEY))
			attrs.put(GaeaVCFConstants.AS_QUAL_KEY, newVC.getAttribute(GaeaVCFConstants.AS_QUAL_KEY));

		return new VariantContextBuilder(newVC).attributes(attrs).make();
	}

	private List<Genotype> cleanupGenotypeAnnotations(final VariantContext VC, final boolean createRefGTs) {
		final GenotypesContext oldGTs = VC.getGenotypes();
		final List<Genotype> recoveredGs = new ArrayList<>(oldGTs.size());
		for (final Genotype oldGT : oldGTs) {
			final Map<String, Object> attrs = new HashMap<>(oldGT.getExtendedAttributes());

			final GenotypeBuilder builder = new GenotypeBuilder(oldGT);
			int depth = oldGT.hasDP() ? oldGT.getDP() : 0;

			// move the MIN_DP to DP
			if (oldGT.hasExtendedAttribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY)) {
				depth = Integer.parseInt((String) oldGT.getAnyAttribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY));
				builder.DP(depth);
				attrs.remove(GaeaVCFConstants.MIN_DP_FORMAT_KEY);
			}

			// move the GQ to RGQ
			if (createRefGTs && oldGT.hasGQ()) {
				builder.noGQ();
				attrs.put(GaeaVCFConstants.REFERENCE_GENOTYPE_QUALITY, oldGT.getGQ());
			}

			// remove SB
			attrs.remove(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);

			// update PGT for hom vars
			if (oldGT.isHomVar() && oldGT.hasExtendedAttribute(GaeaVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY)) {
				attrs.put(GaeaVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "1|1");
			}

			// create AD if it's not there
			if (!oldGT.hasAD() && VC.isVariant()) {
				final int[] AD = new int[VC.getNAlleles()];
				AD[0] = depth;
				builder.AD(AD);
			}

			if (createRefGTs) {
				final int ploidy = oldGT.getPloidy();
				final List<Allele> refAlleles = Collections.nCopies(ploidy, VC.getReference());

				// keep 0 depth samples and 0 GQ samples as no-call
				if (depth > 0 && oldGT.hasGQ() && oldGT.getGQ() > 0) {
					builder.alleles(refAlleles);
				}

				// also, the PLs are technically no longer usable
				builder.noPL();
			}

			recoveredGs.add(builder.noAttributes().attributes(attrs).make());
		}
		return recoveredGs;
	}

	public VCFHeader getVCFHeader() {
		return this.vcfHeader;
	}
}
