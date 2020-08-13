package org.bgi.flexlab.gaea.tools.mapreduce.combinegvcfs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.JointCallingEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.UnifiedGenotypingEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.VariantAnnotatorEngine;
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
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.GaeaVariantContextUtils;
import org.seqdoop.hadoop_bam.LazyParsingGenotypesContext;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.seqdoop.hadoop_bam.LazyVCFGenotypesContext.HeaderDataCache;

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


public class CombineGVCFsEngine extends JointCallingEngine {
	private static String GVCF_BLOCK = "GVCFBlock";
	private final List<String> infoFieldAnnotationKeyNamesToRemove = new ArrayList<>();
	
	private boolean INCLUDE_NON_VARIANTS = false;

	private boolean uniquifySamples = false;

	// private ArrayList<VariantContext> variants = null;
	private TreeMap<String, ArrayList<VariantContext>> variantsForSample = null;
	private String[] samples = null;

	private VariantContext currentContext = null;

	private int max_position = -1;

	// the genotyping engine
	private UnifiedGenotypingEngine genotypingEngine;
	// the annotation engine
	private VariantAnnotatorEngine annotationEngine;

	private GenomeLocationParser parser = null;
	private final VCFHeader vcfHeader;
	protected List<String> annotationsToUse = new ArrayList<>();

	protected List<String> annotationGroupsToUse = new ArrayList<>(
			Arrays.asList(new String[] { StandardAnnotation.class.getSimpleName() }));

	//private final VCFHeader vcfHeader;
	private HashMap<String,HeaderDataCache> vcfHeaderDateCaches = new HashMap<String,HeaderDataCache>();
	
	final Set<String> infoHeaderAltAllelesLineNames = new LinkedHashSet<>();
	public List<VariantContext> getValues(GenomeLocation loc) {
		List<VariantContext> list = new ArrayList<VariantContext>();
		
		if(this.samples != null){
			for (String sample : samples) {
				VariantContext vc = getValues(sample,loc);
				if(vc != null)
					list.add(vc);
			}
		}else{
			for (String sample : variantsForSample.keySet()) {
				VariantContext vc = getValues(sample,loc);
				if(vc != null)
					list.add(vc);
			}
		}
		return list;
	}
	private VariantContext getValues(String sample,GenomeLocation loc,boolean requireStartHere){
		if(variantsForSample.containsKey(sample) && variantsForSample.get(sample).size() > 0){
			for(VariantContext vc : variantsForSample.get(sample)){
				if ( ! requireStartHere || vc.getStart() == loc.getStart()){
					return vc;
				}
			}
		}
		
		return null;
	}
	private VariantContext getValues(String sample,GenomeLocation loc){
		VariantContext vc = getValues(sample,loc,true);
		
		if(vc == null)
			vc = getValues(sample,loc,false);
		
		return vc;
	}
	public CombineGVCFsEngine(JointCallingOptions options, GenomeLocationParser parser, VCFHeader vcfheader,
			MultipleVCFHeaderForJointCalling multiHeaders, String[] sampleArray) {
		
		super(options, parser, multiHeaders);
		variantsForSample = new TreeMap<String, ArrayList<VariantContext>>();
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
		
		headerLines.addAll(vcfheader.getMetaDataInInputOrder());

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

//		for(String sample : multiHeaders.keySet()){
//			HeaderDataCache vcfHeaderDataCache = new HeaderDataCache();
//			vcfHeaderDataCache.setHeader(multiHeaders.getVCFHeader(sample));
//			vcfHeaderDateCaches.put(sample, vcfHeaderDataCache);
//		}

		// now that we have all the VCF headers, initialize the annotations
		// (this is particularly important to turn off RankSumTest dithering in
		// integration tests)
		Set<String> sampleNamesHashSet = new HashSet<>();
		sampleNamesHashSet.addAll(Arrays.asList(sampleArray));
		annotationEngine.invokeAnnotationInitializationMethods(headerLines, sampleNamesHashSet);

		GvcfMathUtils.resetRandomGenerator();
		
		this.samples = sampleArray;
	}
	
	private boolean containsTrueAltAllele(final List<VariantContext> VCs) {
        if ( VCs == null ) throw new IllegalArgumentException("The list of VariantContexts cannot be null");

        for ( final VariantContext vc : VCs ) {
            if ( vc.getNAlleles() > 2 )
                return true;
        }
        return false;
    }
	
	public VariantContext variantCalling3(Iterator<VariantContext> iterator, GenomeLocation location,
			ChromosomeInformationShare ref) {
		
		//lazyLoad(iterator, location);
		if(!iterator.hasNext()) {
			return null;
		}
		//lazyLoad2(iterator,location);
		final List<VariantContext> stoppedVCs = new ArrayList<>();
		int region_start=location.getStart();
		int region_end=location.getEnd();
		if(region_end<region_start) {
			System.exit(1);
		}
		VariantContext tmp_vc;
		while(iterator.hasNext()) {
			//System.out.println("has next:\t"+iterator.next());
			
			tmp_vc=iterator.next();
			String sampleName = tmp_vc.getAttributeAsString("SM", null);
			if (sampleName == null)
				throw new RuntimeException("Not contains SM attribute");
			GenotypesContext gc = tmp_vc.getGenotypes();
			if (gc instanceof LazyParsingGenotypesContext)
				((LazyParsingGenotypesContext) gc).getParser().setHeaderDataCache(vcfHeaderDateCaches.get(sampleName));
			if(tmp_vc.getStart()<=region_end && tmp_vc.getEnd()>=region_start) {
				stoppedVCs.add(tmp_vc);
			}else {
				break;
			}
		}
		if(region_start<20000)
			System.out.println("region start\t"+region_start+"\t"+region_end+"\t"+stoppedVCs.size());
		final Byte refBase = INCLUDE_NON_VARIANTS ? (byte) ref.getBase(location.getStart() - 1) : null;
		final boolean removeNonRefSymbolicAllele = !INCLUDE_NON_VARIANTS;
		VariantContext combinedVC = null;
		if(!stoppedVCs.isEmpty()) {
			if (containsTrueAltAllele(stoppedVCs)) {
				// combinedVC = ReferenceConfidenceVariantContextMerger.merge(stoppedVCs, gLoc,
				// refBase, false, false, annotationEngine);
				System.out.println("variants:\t"+stoppedVCs.toString());
				combinedVC = ReferenceConfidenceVariantContextMerger.merge(stoppedVCs, location, refBase,
						removeNonRefSymbolicAllele, uniquifySamples, annotationEngine);
			} else {
				System.out.println("non variants:\t"+stoppedVCs.toString());
				combinedVC = referenceBlockMerge(stoppedVCs, location);
			}
		}
			
		return combinedVC == null ? null : regenotypeVC(new RefMetaDataTracker(location), ref, combinedVC);
	}
	protected boolean USE_BP_RESOLUTION = false;
	private VariantContext referenceBlockMerge(List<VariantContext> VCs, GenomeLocation location) {
		// TODO Auto-generated method stub
		final VariantContext first = VCs.get(0);
		final Allele refAllele;
        final int start;
        start=location.getStart();
        final int end=location.getEnd();
        refAllele = first.getReference();
        final Map<String, Object> attrs = new HashMap<>(1);
        if ( !USE_BP_RESOLUTION && end != start )
            attrs.put(VCFConstants.END_KEY, Integer.toString(end));
        final GenotypesContext genotypes = GenotypesContext.create();
        for ( final VariantContext vc : VCs ) {
            for ( final Genotype g : vc.getGenotypes() )
                genotypes.add(new GenotypeBuilder(g).alleles(GaeaGvcfVariantContextUtils.noCallAlleles(g.getPloidy())).make());
        }
        return new VariantContextBuilder("", first.getChr(), start, end, Arrays.asList(refAllele, GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE)).attributes(attrs).genotypes(genotypes).make();
	}
	public void lazyLoad2(Iterator<VariantContext> iterator, GenomeLocation location) {
		int curr = location.getStart();

		if (curr <= max_position)
			purgeOutOfScopeRecords(location);
		else {
			for (String sample : variantsForSample.keySet()) {
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
//			if (currentContext.getStart() > curr)
//				break;
			 {
				GenotypesContext gc = currentContext.getGenotypes();
				String sampleName = currentContext.getAttributeAsString("SM", null);
				if (sampleName == null)
					throw new RuntimeException("Not contains SM attribute");
				
				if (gc instanceof LazyParsingGenotypesContext)
					((LazyParsingGenotypesContext) gc).getParser().setHeaderDataCache(vcfHeaderDateCaches.get(sampleName));
				
				if (variantsForSample.containsKey(sampleName)) {
					variantsForSample.get(sampleName).add(currentContext);
				} else {
					ArrayList<VariantContext> list = new ArrayList<VariantContext>();
					list.add(currentContext);
					variantsForSample.put(sampleName, list);
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
	
}
