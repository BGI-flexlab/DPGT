package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.text.SimpleDateFormat;
import java.util.*;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.ReducibleAnnotationData;
import org.bgi.flexlab.gaea.tools.jointcalling.VariantAnnotatorEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculator;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculators;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.GaeaVariantContextUtils;
import org.bgi.flexlab.gaea.util.Pair;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;

public class ReferenceConfidenceVariantContextMerger {

	public static <T extends Comparable<? super T>> T median(final List<T> array) {
		/*
		 * TODO the current implementation is not the usual median when the
		 * input is of even length. More concretely it returns the ith element
		 * of the list where i = floor(input.size() / 2).
		 * 
		 * But actually that is not the "usual" definition of a median, as it is
		 * supposed to return the average of the two middle values when the
		 * sample length is an even number (i.e. median(1,2,3,4,5,6) == 3.5).
		 * [Sources: R and wikipedia]
		 * 
		 * suggestion for a solution is then:
		 * 
		 * unify median and medianDoubles to public static <T extends Number> T
		 * median(Collection<T>) check on null elements and throw an exception
		 * if there are any or perhaps return a null; documented in the javadoc.
		 * relocate, rename and refactor MathUtils.median(X) to
		 * Utils.ithElement(X,X.size()/2) In addition, the current median
		 * implementation sorts the whole input list witch is O(n log n).
		 * However find out the ith element (thus calculate the median) can be
		 * done in O(n)
		 */
		if (array == null)
			throw new IllegalArgumentException("Array must be non-null");
		final int size = array.size();
		if (size == 0)
			throw new IllegalArgumentException("Array cannot have size 0");
		else if (size == 1)
			return array.get(0);
		else {
			final ArrayList<T> sorted = new ArrayList<>(array);
			Collections.sort(sorted);
			return sorted.get(size / 2);
		}
	}

	@SuppressWarnings("unchecked")
	private static Comparable combineAnnotationValues(final List<Comparable> array) {
		return median(array); // right now we take the median but other options
								// could be explored
	}

	public static VariantContext merge(final List<VariantContext> VCs, final GenomeLocation loc, final Byte refBase,
			final boolean removeNonRefSymbolicAllele, final boolean samplesAreUniquified,
			final VariantAnnotatorEngine annotatorEngine) {
		// this can happen if e.g. you are using a dbSNP file that spans a
		// region with no gVCFs
		if (VCs == null || VCs.isEmpty()) {
			return null;
		}
		SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
		// establish the baseline info (sometimes from the first VC)
		final VariantContext first = VCs.get(0);
		final String name = first.getSource();

		// ref allele
		final Allele refAllele = determineReferenceAlleleGivenReferenceBase(VCs, loc, refBase);
		if (refAllele == null) {
			return null;
		}

		// FinalAlleleSet contains the alleles of the new resulting VC
		// Using linked set in order to guarantee a stable order
		final LinkedHashSet<Allele> finalAlleleSet = new LinkedHashSet<>(10);
		// Reference goes first
		finalAlleleSet.add(refAllele);

		final Map<String, Object> attributes = new LinkedHashMap<>();
		final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time
															// there's one id
		int depth = 0;
		final Map<String, List<ReducibleAnnotationData>> annotationMap = new LinkedHashMap<>();
		final GenotypesContext genotypes = GenotypesContext.create();

		// In this list we hold the mapping of each variant context alleles.
		final List<Pair<VariantContext, List<Allele>>> vcAndNewAllelePairs = new ArrayList<>(VCs.size());
		// Keep track of whether we saw a spanning deletion and a non-spanning
		// event
		boolean sawSpanningDeletion = false;
		boolean sawNonSpanningEvent = false;
		// cycle through and add info from the other VCs
		for (final VariantContext vc : VCs) {

			// if this context doesn't start at the current location then it
			// must be a spanning event (deletion or ref block)
			final boolean isSpanningEvent = loc.getStart() != vc.getStart();
			// record whether it's also a spanning deletion/event (we know this
			// because the VariantContext type is no
			// longer "symbolic" but "mixed" because there are real alleles
			// mixed in with the symbolic non-ref allele)
			boolean s = (isSpanningEvent && vc.isMixed())
					|| vc.getAlternateAlleles().contains(Allele.SPAN_DEL)
					|| vc.getAlternateAlleles().contains(GaeaVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
			sawSpanningDeletion |= (isSpanningEvent && vc.isMixed())
					|| vc.getAlternateAlleles().contains(Allele.SPAN_DEL)
					|| vc.getAlternateAlleles().contains(GaeaVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
			sawNonSpanningEvent |= (!isSpanningEvent && vc.isMixed());

			vcAndNewAllelePairs.add(new Pair<>(vc,
					isSpanningEvent ? replaceWithNoCallsAndDels(vc) : remapAlleles(vc, refAllele, finalAlleleSet)));
		}
		// Add <DEL> and <NON_REF> to the end if at all required in the output.
		if (sawSpanningDeletion && (sawNonSpanningEvent || !removeNonRefSymbolicAllele)) {
			finalAlleleSet.add(Allele.SPAN_DEL);
		}
		if (!removeNonRefSymbolicAllele)
			finalAlleleSet.add(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE);

		final List<Allele> allelesList = new ArrayList<>(finalAlleleSet);

		boolean shouldComputePLs = allelesList
				.size() <= GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;
		for (final Pair<VariantContext, List<Allele>> pair : vcAndNewAllelePairs) {
			final VariantContext vc = pair.getFirst();
			final List<Allele> remappedAlleles = pair.getSecond();

			mergeRefConfidenceGenotypes(genotypes, vc, remappedAlleles, allelesList, samplesAreUniquified,
					shouldComputePLs);

			// special case DP (add it up) for all events
			if (vc.hasAttribute(VCFConstants.DEPTH_KEY)) {
				depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
			} else { // handle the gVCF case from the HaplotypeCaller
				for (final Genotype gt : vc.getGenotypes()) {
					depth += (gt.hasExtendedAttribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY)
							? Integer.parseInt((String) gt.getAnyAttribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY))
							: (gt.hasDP() ? gt.getDP() : 0));
				}
			}

			if (loc.getStart() != vc.getStart()) {
				continue;
			}

			// special case ID (just preserve it)
			if (vc.hasID())
				rsIDs.add(vc.getID());

			// add attributes to annotationMap, store all info field annotations
			// as AlleleSpecificAnnotationData in case they can be parsed that
			// way
			addReferenceConfidenceAttributes(pair, annotationMap);
		}
		// combine the annotations that are reducible and remove them from
		// annotationMap
		Map<String, Object> combinedAnnotations = new HashMap<>();
		if (annotatorEngine != null) {
			combinedAnnotations = annotatorEngine.combineAnnotations(allelesList, annotationMap);
		}
		attributes.putAll(combinedAnnotations);

		// remove stale AC and AF based attributes (including MLEAC and MLEAF
		// lists)
		// these will be recalculated after genotyping
		removeStaleAttributesAfterMerge(annotationMap);

		// annotatorEngine.combineAnnotations removed the successfully combined
		// annotations, so now parse those that are left
		// here we're assuming that things that are left are scalars per sample
		Map<String, List<Comparable>> parsedAnnotationMap = parseRemainingAnnotations(annotationMap);

		// when combining remaining annotations use the median value from all
		// input VCs which had annotations provided
		for (final Map.Entry<String, List<Comparable>> p : parsedAnnotationMap.entrySet()) {
			if (!p.getValue().isEmpty()) {
				//System.out.println(p.getKey()+p.getValue().toString());
				attributes.put(p.getKey(), combineAnnotationValues(p.getValue()));
			}
		}

		if (depth > 0) {
			attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));
		}

		final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

		// note that in order to calculate the end position, we need a list of
		// alleles that doesn't include anything symbolic
		final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID).alleles(allelesList)
				.chr(loc.getContig()).start(loc.getStart())
				.computeEndFromAlleles(nonSymbolicAlleles(allelesList), loc.getStart(), loc.getStart())
				.genotypes(genotypes).unfiltered().attributes(new TreeMap<>(attributes))
				.log10PError(CommonInfo.NO_LOG10_PERROR); // we will need to
															// re-genotype later
		return builder.make();
	}
	public static VariantContext mapMerge(final List<VariantContext> VCs, final GenomeLocation loc, final Byte refBase,
									   final boolean removeNonRefSymbolicAllele, final boolean samplesAreUniquified,
									   final VariantAnnotatorEngine annotatorEngine) {
		// this can happen if e.g. you are using a dbSNP file that spans a
		// region with no gVCFs
		if (VCs == null || VCs.isEmpty()) {
			return null;
		}
		SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
		// establish the baseline info (sometimes from the first VC)
		final VariantContext first = VCs.get(0);
		final String name = first.getSource();

		// ref allele
		final Allele refAllele = determineReferenceAlleleGivenReferenceBase(VCs, loc, refBase);
		if (refAllele == null) {
			return null;
		}

		// FinalAlleleSet contains the alleles of the new resulting VC
		// Using linked set in order to guarantee a stable order
		final LinkedHashSet<Allele> finalAlleleSet = new LinkedHashSet<>(10);
		// Reference goes first
		finalAlleleSet.add(refAllele);

		final Map<String, Object> attributes = new LinkedHashMap<>();
		final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time
		// there's one id
		int depth = 0;
		final Map<String, List<ReducibleAnnotationData>> annotationMap = new LinkedHashMap<>();
		final GenotypesContext genotypes = GenotypesContext.create();

		// In this list we hold the mapping of each variant context alleles.
		final List<Pair<VariantContext, List<Allele>>> vcAndNewAllelePairs = new ArrayList<>(VCs.size());
		// Keep track of whether we saw a spanning deletion and a non-spanning
		// event
		boolean sawSpanningDeletion = false;
		boolean sawNonSpanningEvent = false;
		// cycle through and add info from the other VCs
		for (final VariantContext vc : VCs) {

			// if this context doesn't start at the current location then it
			// must be a spanning event (deletion or ref block)
			final boolean isSpanningEvent = loc.getStart() != vc.getStart();
			// record whether it's also a spanning deletion/event (we know this
			// because the VariantContext type is no
			// longer "symbolic" but "mixed" because there are real alleles
			// mixed in with the symbolic non-ref allele)
			boolean s = (isSpanningEvent && vc.isMixed())
					|| vc.getAlternateAlleles().contains(Allele.SPAN_DEL)
					|| vc.getAlternateAlleles().contains(GaeaVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
			sawSpanningDeletion |= (isSpanningEvent && vc.isMixed())
					|| vc.getAlternateAlleles().contains(Allele.SPAN_DEL)
					|| vc.getAlternateAlleles().contains(GaeaVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
			sawNonSpanningEvent |= (!isSpanningEvent && vc.isMixed());

			vcAndNewAllelePairs.add(new Pair<>(vc,
					isSpanningEvent ? replaceWithNoCallsAndDels(vc) : remapAlleles(vc, refAllele, finalAlleleSet)));
		}
		// Add <DEL> and <NON_REF> to the end if at all required in the output.
		if (sawSpanningDeletion && (sawNonSpanningEvent || !removeNonRefSymbolicAllele)) {
			finalAlleleSet.add(Allele.SPAN_DEL);
		}
		if (!removeNonRefSymbolicAllele)
			finalAlleleSet.add(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE);

		final List<Allele> allelesList = new ArrayList<>(finalAlleleSet);

		boolean shouldComputePLs = allelesList
				.size() <= GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;
		HashMap<String,Integer> appearedSamples=new HashMap<String,Integer>();
		for (final Pair<VariantContext, List<Allele>> pair : vcAndNewAllelePairs) {
			final VariantContext vc = pair.getFirst();
			final List<Allele> remappedAlleles = pair.getSecond();

			mergeRefConfidenceGenotypes(genotypes, vc, remappedAlleles, allelesList, samplesAreUniquified,
					shouldComputePLs);
			String sampleName=vc.getAttributeAsString("SM","");
//			if(vc.getStart()==1002417){
//				System.out.println(sampleName);
//				for(String sm:appearedSamples.keySet()){
//					System.out.println(sm+"\t"+appearedSamples.get(sm));
//				}
//			}
			if(appearedSamples.containsKey(sampleName)){
				depth-=appearedSamples.get(sampleName);
				if(depth<0){
					System.out.println("code error\n");
					System.exit(-1);
				}
				appearedSamples.remove(sampleName);
			}
			// special case DP (add it up) for all events
			if (vc.hasAttribute(VCFConstants.DEPTH_KEY)) {
				depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
				appearedSamples.put(sampleName,vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));
			} else { // handle the gVCF case from the HaplotypeCaller
				for (final Genotype gt : vc.getGenotypes()) {
					depth += (gt.hasExtendedAttribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY)
							? Integer.parseInt((String) gt.getAnyAttribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY))
							: (gt.hasDP() ? gt.getDP() : 0));
					if(gt.hasExtendedAttribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY)){
						appearedSamples.put(sampleName,Integer.parseInt((String) gt.getAnyAttribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY)));
					}else{
						if(gt.hasDP()){
							appearedSamples.put(sampleName,gt.getDP());
						}
					}
				}
			}
//			if(loc.getStart()==1002417){
//				System.out.println(vc+"\n"+depth);
//			}
			if (loc.getStart() != vc.getStart()) {
				continue;
			}

			// special case ID (just preserve it)
			if (vc.hasID())
				rsIDs.add(vc.getID());

			// add attributes to annotationMap, store all info field annotations
			// as AlleleSpecificAnnotationData in case they can be parsed that
			// way
			addReferenceConfidenceAttributes(pair, annotationMap);
		}
		// combine the annotations that are reducible and remove them from
		// annotationMap
		Map<String, Object> combinedAnnotations = new HashMap<>();
		if (annotatorEngine != null) {
			combinedAnnotations = annotatorEngine.combineAnnotations(allelesList, annotationMap);
		}
		attributes.putAll(combinedAnnotations);

		// remove stale AC and AF based attributes (including MLEAC and MLEAF
		// lists)
		// these will be recalculated after genotyping
		removeStaleAttributesAfterMerge(annotationMap);

		// annotatorEngine.combineAnnotations removed the successfully combined
		// annotations, so now parse those that are left
		// here we're assuming that things that are left are scalars per sample
		Map<String, List<Comparable>> parsedAnnotationMap = parseRemainingAnnotations(annotationMap);

		// when combining remaining annotations use the median value from all
		// input VCs which had annotations provided
		for (final Map.Entry<String, List<Comparable>> p : parsedAnnotationMap.entrySet()) {
			if (!p.getValue().isEmpty()) {
				//System.out.println(p.getKey()+p.getValue().toString());
				//attributes.put(p.getKey(), combineAnnotationValues(p.getValue()));
				attributes.put(p.getKey(), p.getValue());
			}
		}
		//分样本combine时，可能会出现dp为0的变异，所以这里加上等号
		if (depth >=0) {
			attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));
		}

		final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

		// note that in order to calculate the end position, we need a list of
		// alleles that doesn't include anything symbolic
		final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID).alleles(allelesList)
				.chr(loc.getContig()).start(loc.getStart())
				.computeEndFromAlleles(nonSymbolicAlleles(allelesList), loc.getStart(), loc.getStart())
				.genotypes(genotypes).unfiltered().attributes(new TreeMap<>(attributes))
				.log10PError(CommonInfo.NO_LOG10_PERROR); // we will need to
		// re-genotype later
		return builder.make();
	}
	public static VariantContext reduceMerge(final List<VariantContext> VCs, final GenomeLocation loc, final Byte refBase,
									   final boolean removeNonRefSymbolicAllele, final boolean samplesAreUniquified,
									   final VariantAnnotatorEngine annotatorEngine) {
		// this can happen if e.g. you are using a dbSNP file that spans a
		// region with no gVCFs
		if (VCs == null || VCs.isEmpty()) {
			return null;
		}
		SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
		// establish the baseline info (sometimes from the first VC)
		final VariantContext first = VCs.get(0);
		final String name = first.getSource();

		// ref allele
		final Allele refAllele = determineReferenceAlleleGivenReferenceBase(VCs, loc, refBase);
		if (refAllele == null) {
			return null;
		}

		// FinalAlleleSet contains the alleles of the new resulting VC
		// Using linked set in order to guarantee a stable order
		final LinkedHashSet<Allele> finalAlleleSet = new LinkedHashSet<>(10);
		// Reference goes first
		finalAlleleSet.add(refAllele);

		final Map<String, Object> attributes = new LinkedHashMap<>();
		final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time
		// there's one id
		int depth = 0;
		final Map<String, List<ReducibleAnnotationData>> annotationMap = new LinkedHashMap<>();
		final GenotypesContext genotypes = GenotypesContext.create();

		// In this list we hold the mapping of each variant context alleles.
		final List<Pair<VariantContext, List<Allele>>> vcAndNewAllelePairs = new ArrayList<>(VCs.size());
		// Keep track of whether we saw a spanning deletion and a non-spanning
		// event
		boolean sawSpanningDeletion = false;
		boolean sawNonSpanningEvent = false;
		// cycle through and add info from the other VCs
		for (final VariantContext vc : VCs) {

			// if this context doesn't start at the current location then it
			// must be a spanning event (deletion or ref block)
			final boolean isSpanningEvent = loc.getStart() != vc.getStart();
			// record whether it's also a spanning deletion/event (we know this
			// because the VariantContext type is no
			// longer "symbolic" but "mixed" because there are real alleles
			// mixed in with the symbolic non-ref allele)
			boolean s = (isSpanningEvent && vc.isMixed())
					|| vc.getAlternateAlleles().contains(Allele.SPAN_DEL)
					|| vc.getAlternateAlleles().contains(GaeaVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
			sawSpanningDeletion |= (isSpanningEvent && vc.isMixed())
					|| vc.getAlternateAlleles().contains(Allele.SPAN_DEL)
					|| vc.getAlternateAlleles().contains(GaeaVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
			sawNonSpanningEvent |= (!isSpanningEvent && vc.isMixed());

			vcAndNewAllelePairs.add(new Pair<>(vc,
					isSpanningEvent ? replaceWithNoCallsAndDels(vc) : remapAlleles(vc, refAllele, finalAlleleSet)));
		}
		// Add <DEL> and <NON_REF> to the end if at all required in the output.
		if (sawSpanningDeletion && (sawNonSpanningEvent || !removeNonRefSymbolicAllele)) {
			finalAlleleSet.add(Allele.SPAN_DEL);
		}
		if (!removeNonRefSymbolicAllele)
			finalAlleleSet.add(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE);

		final List<Allele> allelesList = new ArrayList<>(finalAlleleSet);

		boolean shouldComputePLs = allelesList
				.size() <= GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;

		for (final Pair<VariantContext, List<Allele>> pair : vcAndNewAllelePairs) {
			final VariantContext vc = pair.getFirst();
			final List<Allele> remappedAlleles = pair.getSecond();
			mergeRefConfidenceGenotypes(genotypes, vc, remappedAlleles, allelesList, samplesAreUniquified,
					shouldComputePLs);

			// special case DP (add it up) for all events

			if (vc.hasAttribute(VCFConstants.DEPTH_KEY)) {
				if(vc.getAttribute(VCFConstants.DEPTH_KEY,0) instanceof ArrayList){
					System.out.println(vc.getAttribute(VCFConstants.DEPTH_KEY,0));
				}else {
					if(vc.getAttribute(VCFConstants.DEPTH_KEY,0) instanceof String) {
						depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
					}else{
						System.out.println(vc.getAttribute(VCFConstants.DEPTH_KEY,0).getClass());
					}
				}
			} else { // handle the gVCF case from the HaplotypeCaller
				for (final Genotype gt : vc.getGenotypes()) {
					depth += (gt.hasExtendedAttribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY)
							? Integer.parseInt((String) gt.getAnyAttribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY))
							: (gt.hasDP() ? gt.getDP() : 0));
				}
			}

			if (loc.getStart() != vc.getStart()) {
				continue;
			}

			// special case ID (just preserve it)
			if (vc.hasID())
				rsIDs.add(vc.getID());

			// add attributes to annotationMap, store all info field annotations
			// as AlleleSpecificAnnotationData in case they can be parsed that
			// way
			addReferenceConfidenceAttributes(pair, annotationMap);
		}
		// combine the annotations that are reducible and remove them from
		// annotationMap
		Map<String, Object> combinedAnnotations = new HashMap<>();
		if (annotatorEngine != null) {
			combinedAnnotations = annotatorEngine.combineAnnotations(allelesList, annotationMap);
		}
		attributes.putAll(combinedAnnotations);

		// remove stale AC and AF based attributes (including MLEAC and MLEAF
		// lists)
		// these will be recalculated after genotyping
		removeStaleAttributesAfterMerge(annotationMap);

		// annotatorEngine.combineAnnotations removed the successfully combined
		// annotations, so now parse those that are left
		// here we're assuming that things that are left are scalars per sample
		Map<String, List<Comparable>> parsedAnnotationMap = parseRemainingAnnotations2(annotationMap);

		// when combining remaining annotations use the median value from all
		// input VCs which had annotations provided
		for (final Map.Entry<String, List<Comparable>> p : parsedAnnotationMap.entrySet()) {
			if (!p.getValue().isEmpty()) {
				//System.out.println(p.getKey()+p.getValue().toString());
				attributes.put(p.getKey(), combineAnnotationValues(p.getValue()));
			}
		}

		if (depth > 0) {
			attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));
		}

		final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

		// note that in order to calculate the end position, we need a list of
		// alleles that doesn't include anything symbolic
		final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID).alleles(allelesList)
				.chr(loc.getContig()).start(loc.getStart())
				.computeEndFromAlleles(nonSymbolicAlleles(allelesList), loc.getStart(), loc.getStart())
				.genotypes(genotypes).unfiltered().attributes(new TreeMap<>(attributes))
				.log10PError(CommonInfo.NO_LOG10_PERROR); // we will need to
		// re-genotype later
		return builder.make();
	}
	/**
	 * parse the annotations that were not identified as reducible annotations
	 * and combined by the annotation engine
	 * 
	 * @param annotationMap
	 *            the map of info field annotation names and the list of their
	 *            data from the merged VCs
	 * @return info field data parsed as ints or doubles
	 */
	private static Map<String, List<Comparable>> parseRemainingAnnotations(
			final Map<String, List<ReducibleAnnotationData>> annotationMap) {
		final Map<String, List<Comparable>> parsedAnnotations = new HashMap<>();
		for (Map.Entry<String, List<ReducibleAnnotationData>> currentData : annotationMap.entrySet()) {
			List<Comparable> annotationValues = new ArrayList<>();
			if(currentData.getKey().equals("SM")){
				continue;
			}
			for (ReducibleAnnotationData value : currentData.getValue()) {
				try {
					final String stringValue = value.getRawData();
					if (stringValue.contains(".")) {
						annotationValues.add(Double.parseDouble(stringValue));
					} else if (Character.isDigit(stringValue.charAt(0))) {
						if(currentData.getKey().endsWith("Sum")) {
							annotationValues.add(Double.parseDouble(stringValue));
						}else {
							annotationValues.add(Integer.parseInt(stringValue));
						}
						// TODO: uncomment this to parse dbSNP membership
						// annotation once allele-specific merging for that
						// attribute is added
						/*
						 * } else if (Character.isLetter(stringValue.charAt(0)))
						 * { if (stringValue.equalsIgnoreCase("true"))
						 * annotationValues.add(true); else if
						 * (stringValue.equalsIgnoreCase("false"))
						 * annotationValues.add(false);
						 */
					}

				} catch (final NumberFormatException e) {
				}
			}
			parsedAnnotations.put(currentData.getKey(), annotationValues);
		}
		return parsedAnnotations;
	}
	private static Map<String, List<Comparable>> parseRemainingAnnotations2(
			final Map<String, List<ReducibleAnnotationData>> annotationMap) {
		final Map<String, List<Comparable>> parsedAnnotations = new HashMap<>();
		for (Map.Entry<String, List<ReducibleAnnotationData>> currentData : annotationMap.entrySet()) {
			List<Comparable> annotationValues = new ArrayList<>();
			for (ReducibleAnnotationData value : currentData.getValue()) {
				try {
					final String stringValue = value.getRawData();
					if(stringValue.contains(",")){
						String[] eles=stringValue.split(",");
						for(String ele:eles){
							if (ele.contains(".")) {
								annotationValues.add(Double.parseDouble(ele));
							} else if (Character.isDigit(ele.charAt(0))) {
								if (currentData.getKey().endsWith("Sum")) {
									annotationValues.add(Double.parseDouble(ele));
								} else {
									annotationValues.add(Integer.parseInt(ele));
								}
							}
						}
						continue;
					}
					if (stringValue.contains(".")) {
						annotationValues.add(Double.parseDouble(stringValue));
					} else if (Character.isDigit(stringValue.charAt(0))) {
						if(currentData.getKey().endsWith("Sum")) {
							annotationValues.add(Double.parseDouble(stringValue));
						}else {
							annotationValues.add(Integer.parseInt(stringValue));
						}
						// TODO: uncomment this to parse dbSNP membership
						// annotation once allele-specific merging for that
						// attribute is added
						/*
						 * } else if (Character.isLetter(stringValue.charAt(0)))
						 * { if (stringValue.equalsIgnoreCase("true"))
						 * annotationValues.add(true); else if
						 * (stringValue.equalsIgnoreCase("false"))
						 * annotationValues.add(false);
						 */
					}

				} catch (final NumberFormatException e) {
				}
			}
			parsedAnnotations.put(currentData.getKey(), annotationValues);
		}
		return parsedAnnotations;
	}
	/**
	 * Remove the stale attributes from the merged set
	 *
	 * @param attributes
	 *            the attribute map
	 */
	private static void removeStaleAttributesAfterMerge(final Map<String, List<ReducibleAnnotationData>> attributes) {
		attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
		attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
		attributes.remove(VCFConstants.ALLELE_NUMBER_KEY);
		attributes.remove(GaeaVCFConstants.MLE_ALLELE_COUNT_KEY);
		attributes.remove(GaeaVCFConstants.MLE_ALLELE_FREQUENCY_KEY);
		attributes.remove(VCFConstants.END_KEY);
	}

	private static void addReferenceConfidenceAttributes(Pair<VariantContext, List<Allele>> pair,
			final Map<String, List<ReducibleAnnotationData>> annotationMap) {
		final Map<String, Object> myAttributes = pair.getFirst().getAttributes();
		final List<Allele> sampleAlleles = pair.getSecond();

		for (final Map.Entry<String, Object> p : myAttributes.entrySet()) {
			final String key = p.getKey();
			// allele-specific attributes will always be in list form because
			// they've already been parsed per-allele
			// non-allele-specific attributes (DP, etc.) will be a list of
			// length 1
			final List<Object> valueList = pair.getFirst().getAttributeAsList(key);

			// add the existing annotation values to a list for combining later
			List<ReducibleAnnotationData> rawValuesList = annotationMap.get(key);
			if (rawValuesList == null) {
				rawValuesList = new ArrayList<>();
				annotationMap.put(key, rawValuesList);
			}
			String combinedString = "";
			for (int i = 0; i < valueList.size(); i++) {
				if (i > 0)
					combinedString += ",";
				combinedString += valueList.get(i);
			}
			ReducibleAnnotationData pairData = new AlleleSpecificAnnotationData(sampleAlleles, combinedString);
			rawValuesList.add(pairData);
			annotationMap.put(key, rawValuesList);
		}
	}

	/**
	 * @param list
	 *            the original alleles list
	 * @return a non-null list of non-symbolic alleles
	 */
	private static List<Allele> nonSymbolicAlleles(final List<Allele> list) {
		final List<Allele> result = new ArrayList<>(list.size());
		for (final Allele allele : list) {
			if (!allele.isSymbolic()) {
				result.add(allele);
			}
		}
		return result;
	}

	private static Allele determineReferenceAlleleGivenReferenceBase(final List<VariantContext> VCs,
			final GenomeLocation loc, final Byte refBase) {
		final Allele refAllele = GaeaVariantContextUtils.determineReferenceAllele(VCs, loc);
		if (refAllele == null) {
			return (refBase == null ? null : Allele.create(refBase, true));
		}
		return refAllele;
	}

	/**
	 * Replaces any alleles in the VariantContext with NO CALLS or the symbolic
	 * deletion allele as appropriate, except for the generic ALT allele
	 *
	 * @param vc
	 *            VariantContext with the alleles to replace
	 * @return non-null list of alleles
	 */
	private static List<Allele> replaceWithNoCallsAndDels(final VariantContext vc) {
		if (vc == null)
			throw new IllegalArgumentException("VariantContext cannot be null");

		final List<Allele> result = new ArrayList<>(vc.getNAlleles());

		// no-call the reference allele
		result.add(Allele.NO_CALL);

		// handle the alternate alleles
		for (final Allele allele : vc.getAlternateAlleles()) {
			final Allele replacement;
			if (allele.equals(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE))
				replacement = allele;
			else if (allele.length() < vc.getReference().length())
				replacement = Allele.SPAN_DEL;
			else
				replacement = Allele.NO_CALL;

			result.add(replacement);
		}
		return result;
	}

	private static List<Allele> remapAlleles(final VariantContext vc, final Allele refAllele,
			final LinkedHashSet<Allele> finalAlleles) {

		final Allele vcRef = vc.getReference();
		final byte[] refBases = refAllele.getBases();
		final int extraBaseCount = refBases.length - vcRef.getBases().length;
		if (extraBaseCount < 0)
			throw new IllegalStateException("the wrong reference was selected");

		final List<Allele> result = new ArrayList<>(vc.getNAlleles());
		result.add(refAllele);

		for (final Allele a : vc.getAlternateAlleles()) {
			if (a.isSymbolic()) {
				result.add(a);
				// we always skip <NON_REF> when adding to finalAlleles; this is
				// done outside if it applies.
				// we also skip <*:DEL> if there isn't a real alternate allele.
				if (!a.equals(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE) && !vc.isSymbolic())
					finalAlleles.add(a);
			} else if (a == Allele.SPAN_DEL) {
				result.add(a);
				// we skip * if there isn't a real alternate allele.
				if (!vc.isBiallelic())
					finalAlleles.add(a);
			} else if (a.isCalled()) {
				final Allele newAllele;
				if (extraBaseCount > 0) {
					final byte[] oldBases = a.getBases();
					final byte[] newBases = Arrays.copyOf(oldBases, oldBases.length + extraBaseCount);
					System.arraycopy(refBases, refBases.length - extraBaseCount, newBases, oldBases.length,
							extraBaseCount);
//					byte[] tmpBases={'*','T','G'};
//					if(Arrays.equals(newBases,tmpBases)){
//						System.out.println("here");
//					}
					if(newBases[0]=='*'){
						byte[] newBases2=new byte[newBases.length-1];
						for(int i=1;i<newBases.length;i++){
							newBases2[i-1]=newBases[i];
						}
						newAllele = Allele.create(newBases2, false);
					}else {
						newAllele = Allele.create(newBases, false);
					}
				} else
					newAllele = a;
				result.add(newAllele);
				finalAlleles.add(newAllele);
			} else { // NO_CALL and strange miscellanea
				result.add(a);
			}
		}
		return result;
	}

	/**
	 * Merge into the context a new genotype represented by the given
	 * VariantContext for the provided list of target alleles. This method
	 * assumes that none of the alleles in the VC overlaps with any of the
	 * alleles in the set.
	 */
	private static void mergeRefConfidenceGenotypes(final GenotypesContext mergedGenotypes, final VariantContext vc,
			final List<Allele> remappedAlleles, final List<Allele> targetAlleles, final boolean samplesAreUniquified,
			final boolean shouldComputePLs) {
		final int maximumPloidy = vc.getMaxPloidy(GaeaGvcfVariantContextUtils.DEFAULT_PLOIDY);
		// the map is different depending on the ploidy, so in order to keep
		// this method flexible (mixed ploidies)
		// we need to get a map done (lazily inside the loop) for each ploidy,
		// up to the maximum possible.
		final int[][] genotypeIndexMapsByPloidy = new int[maximumPloidy + 1][];
		final int maximumAlleleCount = Math.max(remappedAlleles.size(), targetAlleles.size());

		for (final Genotype g : vc.getGenotypes()) {
			final String name;
			if (samplesAreUniquified)
				name = g.getSampleName() + "." + vc.getSource();
			else
				name = g.getSampleName();
			final int ploidy = g.getPloidy();
			final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(g)
					.alleles(GaeaGvcfVariantContextUtils.noCallAlleles(g.getPloidy())).noPL();
			genotypeBuilder.name(name);

			final boolean doPLs = shouldComputePLs && g.hasPL();
			final boolean hasAD = g.hasAD();
			final boolean hasSAC = g.hasExtendedAttribute(GaeaVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY);
			if (doPLs || hasSAC || hasAD) {
				final int[] perSampleIndexesOfRelevantAlleles = getIndexesOfRelevantAlleles(remappedAlleles,
						targetAlleles, vc.getStart(), g);
				if (doPLs) {
					// lazy initialization of the genotype index map by ploidy.

					final int[] genotypeIndexMapByPloidy = genotypeIndexMapsByPloidy[ploidy] == null
							? GenotypeLikelihoodCalculators.getInstance(ploidy, maximumAlleleCount).genotypeIndexMap(
									perSampleIndexesOfRelevantAlleles)
							: genotypeIndexMapsByPloidy[ploidy];
					final int[] PLs = generatePL(g, genotypeIndexMapByPloidy);
					genotypeBuilder.PL(PLs);
				}
				if (hasAD) {
					genotypeBuilder.AD(generateAD(g.getAD(), perSampleIndexesOfRelevantAlleles));
				}
				if (hasSAC) {
					final List<Integer> sacIndexesToUse = adaptToSACIndexes(perSampleIndexesOfRelevantAlleles);
					final int[] SACs = GaeaGvcfVariantContextUtils.makeNewSACs(g, sacIndexesToUse);
					genotypeBuilder.attribute(GaeaVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, SACs);
				}
			}
			mergedGenotypes.add(genotypeBuilder.make());
		}
	}

	private static int[] generatePL(final Genotype g, final int[] genotypeIndexMapByPloidy) {
		final int[] PLs = new int[genotypeIndexMapByPloidy.length];
		final int[] oldPLs = g.getPL();
		for (int i = 0; i < PLs.length; i++)
			PLs[i] = oldPLs[genotypeIndexMapByPloidy[i]];
		return PLs;
	}

	protected static int[] generateAD(final int[] originalAD, final int[] indexesOfRelevantAlleles) {
		if (originalAD == null || indexesOfRelevantAlleles == null)
			throw new IllegalArgumentException("The list of input AD values and alleles must not be null");

		final int numADs = indexesOfRelevantAlleles.length;
		final int[] newAD = new int[numADs];

		for (int i = 0; i < numADs; i++) {
			final int oldIndex = indexesOfRelevantAlleles[i];
			if (oldIndex >= originalAD.length)
				newAD[i] = 0;
			else
				newAD[i] = originalAD[oldIndex];
		}

		return newAD;
	}

	private static List<Integer> adaptToSACIndexes(final int[] perSampleIndexesOfRelevantAlleles) {
		if (perSampleIndexesOfRelevantAlleles == null)
			throw new IllegalArgumentException("The per sample index of relevant alleles must not be null");

		final List<Integer> sacIndexesToUse = new ArrayList(2 * perSampleIndexesOfRelevantAlleles.length);

		for (int item : perSampleIndexesOfRelevantAlleles) {
			sacIndexesToUse.add(new Integer(2 * item));
			sacIndexesToUse.add(new Integer(2 * item + 1));
		}

		return sacIndexesToUse;
	}

	protected static int[] getIndexesOfRelevantAlleles(final List<Allele> remappedAlleles,
			final List<Allele> targetAlleles, final int position, final Genotype g) {

		if (remappedAlleles == null || remappedAlleles.isEmpty())
			throw new IllegalArgumentException("The list of input alleles must not be null or empty");
		if (targetAlleles == null || targetAlleles.isEmpty())
			throw new IllegalArgumentException("The list of target alleles must not be null or empty");

		if (!remappedAlleles.contains(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE))
			throw new UserException("The list of input alleles must contain " + GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE
					+ " as an allele but that is not the case at position " + position
					+ "; please use the Haplotype Caller with gVCF output to generate appropriate records");

		final int indexOfNonRef = remappedAlleles.indexOf(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE);
		final int[] indexMapping = new int[targetAlleles.size()];

		// the reference likelihoods should always map to each other (even if
		// the alleles don't)
		indexMapping[0] = 0;

		// create the index mapping, using the <NON-REF> allele whenever such a
		// mapping doesn't exist
		for (int i = 1; i < targetAlleles.size(); i++) {
			final Allele targetAllele = targetAlleles.get(i);

			// if there’s more than 1 DEL allele then we need to use the best
			// one
			if (targetAllele == Allele.SPAN_DEL && g.hasPL()) {
				final int occurrences = Collections.frequency(remappedAlleles, Allele.SPAN_DEL);
				if (occurrences > 1) {
					final int indexOfBestDel = indexOfBestDel(remappedAlleles, g.getPL(), g.getPloidy());
					indexMapping[i] = (indexOfBestDel == -1 ? indexOfNonRef : indexOfBestDel);
					continue;
				}
			}

			final int indexOfRemappedAllele = remappedAlleles.indexOf(targetAllele);
			indexMapping[i] = indexOfRemappedAllele == -1 ? indexOfNonRef : indexOfRemappedAllele;
		}

		return indexMapping;
	}

	private static int indexOfBestDel(final List<Allele> alleles, final int[] PLs, final int ploidy) {
		int bestIndex = -1;
		int bestPL = Integer.MAX_VALUE;

		for (int i = 0; i < alleles.size(); i++) {
			if (alleles.get(i) == Allele.SPAN_DEL) {
				final int homAltIndex = findHomIndex(i, ploidy, alleles.size());
				final int PL = PLs[homAltIndex];
				if (PL < bestPL) {
					bestIndex = i;
					bestPL = PL;
				}
			}
		}

		return bestIndex;
	}

	private static int findHomIndex(final int i, final int ploidy, final int numAlleles) {
		// some quick optimizations for the common case
		if (ploidy == 2)
			return i + (i * (i + 1) / 2); // this is straight from the VCF spec
											// on PLs
		if (ploidy == 1)
			return i;

		final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(ploidy, numAlleles);
		final int[] alleleIndexes = new int[ploidy];
		Arrays.fill(alleleIndexes, i);
		return calculator.allelesToIndex(alleleIndexes);
	}
}
