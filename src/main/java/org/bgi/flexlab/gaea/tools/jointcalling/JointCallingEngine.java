package org.bgi.flexlab.gaea.tools.jointcalling;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.input.vcf.VCFHdfsLoader;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.annotator.RMSAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.annotator.RankSumTest;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.bgi.flexlab.gaea.tools.jointcalling.util.ReferenceConfidenceVariantContextMerger;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCallingOptions;
import org.bgi.flexlab.gaea.util.FileIterator;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.seqdoop.hadoop_bam.LazyParsingGenotypesContext;
import org.seqdoop.hadoop_bam.LazyVCFGenotypesContext.HeaderDataCache;
import org.seqdoop.hadoop_bam.VariantContextWithHeader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import org.seqdoop.hadoop_bam.VariantContextWritable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import htsjdk.variant.vcf.VCFUtils;

public class JointCallingEngine {
	
	private static String GVCF_BLOCK = "GVCFBlock";
	private final List<String> infoFieldAnnotationKeyNamesToRemove = new ArrayList<>();
	
	private boolean INCLUDE_NON_VARIANTS = false;

	private boolean uniquifySamples = false;

	// private ArrayList<VariantContext> variants = null;
	private TreeMap<Integer, ArrayList<VariantContext>> variantsForSample = null;
	//private String[] samples = null;

	private VariantContext currentContext = null;
	private static Set<VCFHeaderLine> gvcfHeaderMetaInfo;
	private int max_position = -1;

	// the genotyping engine
	private UnifiedGenotypingEngine genotypingEngine;
	// the annotation engine
	private VariantAnnotatorEngine annotationEngine;

	private GenomeLocationParser parser = null;

	protected List<String> annotationsToUse = new ArrayList<>();

	protected List<String> annotationGroupsToUse = new ArrayList<>(
			Arrays.asList(new String[] { StandardAnnotation.class.getSimpleName() }));
	SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
	private final VCFHeader vcfHeader;
	private HashMap<Integer,HeaderDataCache> vcfHeaderDateCaches = null;
	private HashMap<Integer,Set<String>> nameHeaders=new HashMap<Integer,Set<String>>();
	final Set<String> infoHeaderAltAllelesLineNames = new LinkedHashSet<>();
	private Set<Integer> mapMergedSamples=new TreeSet<Integer>();
	private int sampleSize = 0;
	public JointCallingEngine(JointCallingOptions options, GenomeLocationParser parser,
			MultipleVCFHeaderForJointCalling multiHeaders) {
		variantsForSample = new TreeMap<Integer, ArrayList<VariantContext>>();
		this.INCLUDE_NON_VARIANTS = options.INCLUDE_NON_VARIANT;
		this.uniquifySamples = options.isUniquifySamples();
		this.parser = parser;

		annotationEngine = new VariantAnnotatorEngine(annotationGroupsToUse, annotationsToUse,
				Collections.<String>emptyList());
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
		final Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();

		headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().contains(GVCF_BLOCK));
		
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

		/*for(String sample : multiHeaders.keySet()){
			HeaderLineDataCache vcfHeaderDataCache = new HeaderLineDataCache();
			vcfHeaderDataCache.setHeader(multiHeaders.getVCFHeader(sample));
			vcfHeaderDateCaches.put(sample, vcfHeaderDataCache);
		}*/
		
		this.sampleSize = multiHeaders.getHeaderSize();
		//this.vcfHeaderDateCaches = multiHeaders.getHeaderDataCache();

		// now that we have all the VCF headers, initialize the annotations
		// (this is particularly important to turn off RankSumTest dithering in
		// integration tests)
		Set<String> sampleNamesHashSet = new HashSet<>();
		sampleNamesHashSet.addAll(multiHeaders.getMergeHeader().getSampleNamesInOrder());
		annotationEngine.invokeAnnotationInitializationMethods(headerLines, sampleNamesHashSet);
		gvcfHeaderMetaInfo=multiHeaders.getMergeHeader().getMetaDataInInputOrder();
		GvcfMathUtils.resetRandomGenerator();

	}

	public JointCallingEngine(JointCallingOptions options, GenomeLocationParser parser, VCFHeader vcfheader,
			MultipleVCFHeaderForJointCalling multiHeaders,String[] sampleArray,ArrayList<ArrayList<String> >multiMapSamples,Configuration conf) throws IllegalArgumentException, IOException{
		variantsForSample = new TreeMap<Integer, ArrayList<VariantContext>>();
		this.INCLUDE_NON_VARIANTS = options.INCLUDE_NON_VARIANT;
		this.uniquifySamples = options.isUniquifySamples();
		this.parser = parser;

		annotationEngine = new VariantAnnotatorEngine(annotationGroupsToUse, annotationsToUse,
				Collections.<String>emptyList());
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
		final Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();

		headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().contains(GVCF_BLOCK));

		headerLines.addAll(multiHeaders.getMergeHeader().getMetaDataInInputOrder());


		//headerLines.addAll(vcfheader.getMetaDataInInputOrder());

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
//		
		//this.vcfHeaderDateCaches = multiHeaders.getHeaderDataCache();
		ArrayList<String> totalSamples=new ArrayList<String>();
//		String[] newSamples=new String[multiMapSamples.size()];
//		int i=0;
		HashMap<String, Integer> sampleIndex=new HashMap();
		String sampleIndexFile=options.getOutDir()+"/vcfheaderinfo";
		Path sampleIndexFilePath=new Path(sampleIndexFile);
		FileSystem fs=sampleIndexFilePath.getFileSystem(conf);
		BufferedReader reader=null;
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
			Set<String> samples = new LinkedHashSet<String>();
			String[] sample = eles[2].split(",");
			for (String s : sample) {
				samples.add(s);
			}
			nameHeaders.put(Integer.parseInt(eles[1]), samples);
		}
		reader.close();
		//create nameHeaders
		

		File rawInput=new File(options.getInputList());
		String sortedInputGvcfList=rawInput.getParent()+"/sorted."+rawInput.getName();
		if(sortedInputGvcfList.startsWith("file:")){
			sortedInputGvcfList=sortedInputGvcfList.substring(5);
		}
		BufferedReader inputListReader=null;
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
		String pathLine;
		String[] paths=new String[sampleIndex.size()];
		int i=0;
		while((pathLine=inputListReader.readLine())!=null) {
			paths[i]=pathLine;
			i++;
		}
		inputListReader.close();
		for(ArrayList<String> mapSamples:multiMapSamples) {
			HeaderDataCache vcfHeaderDataCache = new HeaderDataCache();
			VCFHeader MergedMultiSamplesHeader=new VCFHeader();
			Set<VCFHeader> multiSamples=new HashSet<VCFHeader>();
			Set<String> samples = new LinkedHashSet<String>();
			for(String sm:mapSamples) {
				samples.add(sm);
//				String samplePath=paths[sampleIndex.get(sm)];
//				if(samplePath.startsWith("file://")){
//					samplePath=samplePath.substring(7);
//				}
//				VCFLocalLoader curSampleHeaderLoader=new VCFLocalLoader(samplePath);
//				multiSamples.add(curSampleHeaderLoader.getHeader());
//				curSampleHeaderLoader.close();
			}
//			MergedMultiSamplesHeader=smartMergeHeaders(multiSamples);
//			vcfHeaderDataCache.setHeader(MergedMultiSamplesHeader);
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
			//vcfHeaderDateCaches.put(mergedSampleIndex,vcfHeaderDataCache);
		}
		
		// now that we have all the VCF headers, initialize the annotations
		// (this is particularly important to turn off RankSumTest dithering in
		// integration tests)


		Set<String> sampleNamesHashSet = new HashSet<>();
		sampleNamesHashSet.addAll(totalSamples);
		annotationEngine.invokeAnnotationInitializationMethods(headerLines, sampleNamesHashSet);
		gvcfHeaderMetaInfo=multiHeaders.getMergeHeader().getMetaDataInInputOrder();
		GvcfMathUtils.resetRandomGenerator();
//		
//		this.samples=newSamples;
	}
	public static Set<String> getSampleList(Set<VCFHeader> headers) {
		Set<String> samples = new TreeSet<String>();
		for (VCFHeader header : headers) {
			for (String sample : header.getGenotypeSamples()) {
				samples.add(GaeaGvcfVariantContextUtils.mergedSampleName(null, sample, false));
			}
		}

		return samples;
	}
	private VCFHeader smartMergeHeaders(Set<VCFHeader> headers) {
		Set<String> samplelists = getSampleList(headers);
		Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(headers, true);
		VCFHeader vcfHeader = new VCFHeader(headerLines, samplelists);
		headers.clear();
		samplelists.clear();
		headerLines.clear();

		return vcfHeader;
	}
	public Set<String> getSampleList(VCFHeader header) {
		Set<String> samples = new TreeSet<String>();
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
					ArrayList<VariantContext> list = new ArrayList<VariantContext>();
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
					ArrayList<VariantContext> list = new ArrayList<VariantContext>();
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
	public void lazyLoadForSpark2(Iterator<VariantContext> iterator, GenomeLocation location) {
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
					ArrayList<VariantContext> list = new ArrayList<VariantContext>();
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
		List<VariantContext> list = new ArrayList<VariantContext>();
		
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

	public VariantContext variantCalling(Iterator<VariantContextWritable> iterator, GenomeLocation location,
			ChromosomeInformationShare ref) {
		if (location.getStart() != location.getStop())
			throw new UserException("location must length is 1!");
		lazyLoad(iterator, location);
//		for(Entry<Integer, ArrayList<VariantContext> > entry:variantsForSample.entrySet()) {
//			System.out.println(entry.getKey()+"\t"+entry.getValue().size());
//			for(VariantContext tmpVC:entry.getValue()) {
//				System.out.println(tmpVC);
//			}
//		}
		final List<VariantContext> vcsAtThisLocus = getValues(location);
		final Byte refBase = INCLUDE_NON_VARIANTS ? (byte) ref.getBase(location.getStart() - 1) : null;
		final boolean removeNonRefSymbolicAllele = !INCLUDE_NON_VARIANTS;
		final VariantContext combinedVC = ReferenceConfidenceVariantContextMerger.merge(vcsAtThisLocus, location,
				refBase, removeNonRefSymbolicAllele, uniquifySamples, annotationEngine);
		return combinedVC == null ? null : regenotypeVC(new RefMetaDataTracker(location), ref, combinedVC);
	}
	public VariantContext variantCallingReduce(Iterator<VariantContextWritable> iterator, GenomeLocation location,
										 ChromosomeInformationShare ref) {
		if (location.getStart() != location.getStop())
			throw new UserException("location must length is 1!");
		lazyLoad(iterator, location);
//		for(Entry<Integer, ArrayList<VariantContext> > entry:variantsForSample.entrySet()) {
//			System.out.println(entry.getKey()+"\t"+entry.getValue().size());
//			for(VariantContext tmpVC:entry.getValue()) {
//				System.out.println(tmpVC);
//			}
//		}
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
//		for(Entry<Integer, ArrayList<VariantContext> > entry:variantsForSample.entrySet()) {
//			System.out.println(entry.getKey()+"\t"+entry.getValue().size());
//			for(VariantContext tmpVC:entry.getValue()) {
//				System.out.println(tmpVC);
//			}
//		}
		final List<VariantContext> vcsAtThisLocus = getValues(location);
		final Byte refBase = INCLUDE_NON_VARIANTS ? (byte) ref.getBase(location.getStart() - 1) : null;
		final boolean removeNonRefSymbolicAllele = !INCLUDE_NON_VARIANTS;
		final VariantContext combinedVC = ReferenceConfidenceVariantContextMerger.reduceMerge(vcsAtThisLocus, location,
				refBase, removeNonRefSymbolicAllele, uniquifySamples, annotationEngine);
		return combinedVC == null ? null : regenotypeVC(new RefMetaDataTracker(location), ref, combinedVC);
	}
	public VariantContext variantCallingForSpark2(Iterator<VariantContext> iterator, GenomeLocation location,
												 ChromosomeInformationShare ref) {
		if (location.getStart() != location.getStop())
			throw new UserException("location must length is 1!");
		lazyLoadForSpark2(iterator, location);
//		for(Entry<Integer, ArrayList<VariantContext> > entry:variantsForSample.entrySet()) {
//			System.out.println(entry.getKey()+"\t"+entry.getValue().size());
//			for(VariantContext tmpVC:entry.getValue()) {
//				System.out.println(tmpVC);
//			}
//		}
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
                builder.rmAttributes(Arrays.asList(name));
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
