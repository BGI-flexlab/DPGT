package org.bgi.flexlab.gaea.tools.jointcalling;

import java.text.SimpleDateFormat;
import java.util.*;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.variant.VariantCallContext;
import org.bgi.flexlab.gaea.tools.jointcalling.afcalculator.AFCalculationResult;
import org.bgi.flexlab.gaea.tools.jointcalling.afcalculator.AFCalculator;
import org.bgi.flexlab.gaea.tools.jointcalling.afcalculator.AFCalculatorProvider;
import org.bgi.flexlab.gaea.tools.jointcalling.afcalculator.AlleleFrequencyCalculator;
import org.bgi.flexlab.gaea.tools.jointcalling.afcalculator.GeneralPloidyFailOverAFCalculatorProvider;
import org.bgi.flexlab.gaea.tools.jointcalling.afpriorprovider.AFPriorProvider;
import org.bgi.flexlab.gaea.tools.jointcalling.afpriorprovider.CustomAFPriorProvider;
import org.bgi.flexlab.gaea.tools.jointcalling.afpriorprovider.HeterozygosityAFPriorProvider;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCallingOptions;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class UnifiedGenotypingEngine {

	protected final int numberOfGenomes;
	private boolean doAlleleSpecificCalcs = false;
	private GenomeLocationParser genomeLocParser;

	private final AFPriorProvider log10AlleleFrequencyPriorsSNPs;

	private final AFPriorProvider log10AlleleFrequencyPriorsIndels;

	private JointCallingOptions options = null;

	protected final AFCalculatorProvider afCalculatorProvider;

	private final List<GenomeLocation> upstreamDeletionsLoc = new LinkedList<>();

	protected final AFCalculator newAFCalculator;

	public UnifiedGenotypingEngine(int sampleCount, JointCallingOptions options, GenomeLocationParser parser) {
		this.options = options;
		this.genomeLocParser = parser;
		numberOfGenomes = sampleCount * options.getSamplePloidy();
		GvcfMathUtils.Log10Cache.ensureCacheContains(numberOfGenomes * 2);
		log10AlleleFrequencyPriorsSNPs = composeAlleleFrequencyPriorProvider(numberOfGenomes,
				options.getSNPHeterozygosity(), options.getInputPrior());
		log10AlleleFrequencyPriorsIndels = composeAlleleFrequencyPriorProvider(numberOfGenomes,
				options.getINDELHeterozygosity(), options.getInputPrior());

		afCalculatorProvider = new GeneralPloidyFailOverAFCalculatorProvider(options.getSamplePloidy(),
				options.getMaxAlternateAllele());
		
		final double refPseudocount = options.getSNPHeterozygosity() / Math.pow(options.heterozygosityStandardDeviation,2);
	    final double snpPseudocount = options.getSNPHeterozygosity() * refPseudocount;
	    final double indelPseudocount = options.getINDELHeterozygosity() * refPseudocount;
	    newAFCalculator = new AlleleFrequencyCalculator(refPseudocount, snpPseudocount, indelPseudocount, options.getSamplePloidy());
	}

	public enum Model {
		SNP, INDEL, GENERALPLOIDYSNP, GENERALPLOIDYINDEL, BOTH;
	}

	public enum GenotypingOutputMode {

		/**
		 * The genotyper will choose the most likely alternate allele
		 */
		DISCOVERY,

		/**
		 * Only the alleles passed by the user should be considered.
		 */
		GENOTYPE_GIVEN_ALLELES
	}

	public enum OutputMode {
		/** produces calls only at variant sites */
		EMIT_VARIANTS_ONLY,
		/** produces calls at variant sites and confident reference sites */
		EMIT_ALL_CONFIDENT_SITES,
		/**
		 * produces calls at any callable site regardless of confidence; this
		 * argument is intended only for point mutations (SNPs) in DISCOVERY
		 * mode or generally when running in GENOTYPE_GIVEN_ALLELES mode; it
		 * will by no means produce a comprehensive set of indels in DISCOVERY
		 * mode
		 */
		EMIT_ALL_SITES
	}

	private static class OutputAlleleSubset {
		private final Allele[] alleles;
		private final boolean siteIsMonomorphic;
		private final int[] mleCounts;
		private final int count;

		private OutputAlleleSubset(final int count, final Allele[] alleles, final int[] mleCounts,
				final boolean siteIsMonomorphic) {
			Utils.nonNull(alleles, "alleles must not be null");
			Utils.nonNull(mleCounts, "mleCounts must not be null");
			Utils.validateArg(count <= alleles.length, "count must be <= "+alleles.length);
			Utils.validateArg(count <= mleCounts.length, "count must be <= "+mleCounts.length);
			this.siteIsMonomorphic = siteIsMonomorphic;
			this.count = count;
			this.alleles = alleles;
			this.mleCounts = mleCounts;
		}

		private List<Allele> outputAlleles(final Allele referenceAllele) {
			final ArrayList<Allele> result = new ArrayList<>(count + 1);
			result.add(referenceAllele);
			for (int i = 0; i < count; i++)
				result.add(alleles[i]);
			return result;
		}

		public List<Integer> alternativeAlleleMLECounts() {
			final List<Integer> result = new ArrayList<>(count);
			for (int i = 0; i < count; i++)
				result.add(mleCounts[i]);
			return result;
		}
	}

	public static AFPriorProvider composeAlleleFrequencyPriorProvider(final int N, final double heterozygosity,
			final List<Double> inputPriors) {

		if (!inputPriors.isEmpty()) {
			// user-specified priors
			if (inputPriors.size() != N)
				throw new UserException(
						"Invalid length of inputPrior vector: vector length must be equal to # samples +1 ");
			for (final Double prior : inputPriors) {
				if (prior <= 0 || prior >= 1)
					throw new UserException("inputPrior vector values must be greater than 0 and less than 1");
			}
			return new CustomAFPriorProvider(inputPriors);
		} else
			return new HeterozygosityAFPriorProvider(heterozygosity);
	}

	public VariantCallContext calculateGenotypes(VariantContext vc) {
		final VariantContext.Type type = vc.getType();
		final Model model;

		if (type == VariantContext.Type.INDEL) {
			model = Model.INDEL;
		} else {
			model = Model.SNP;
		}

		return calculateGenotypes(vc, model, doAlleleSpecificCalcs);
	}

	public VariantCallContext calculateGenotypes(final VariantContext vc, final Model model,
			final boolean useAlleleSpecificCalcs) {
		return calculateGenotypes(vc, model, false, useAlleleSpecificCalcs);
	}

	protected final boolean hasTooManyAlternativeAlleles(final VariantContext vc) {
		// protect against too many alternate alleles that we can't even run AF
		// on:
		Utils.nonNull(vc, "variantContext must not be null");
		if (vc.getNAlleles() <= GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED)
			return false;
		return true;
	}

	protected boolean forceSiteEmission() {
		return options.getOutputMode() == OutputMode.EMIT_ALL_SITES;
	}

	protected final double[] getAlleleFrequencyPriors(final VariantContext vc, final int defaultPloidy,
			final Model model) {
		final int totalPloidy = GaeaGvcfVariantContextUtils.totalPloidy(vc, defaultPloidy);
		switch (model) {
		case SNP:
		case GENERALPLOIDYSNP:
			return log10AlleleFrequencyPriorsSNPs.forTotalPloidy(totalPloidy);
		case INDEL:
		case GENERALPLOIDYINDEL:
			return log10AlleleFrequencyPriorsIndels.forTotalPloidy(totalPloidy);
		default:
			throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + model);
		}
	}

	protected final boolean confidentlyCalled(final double conf, final double PofF) {
		return conf >= options.STANDARD_CONFIDENCE_FOR_CALLING
				|| (options.getGenotypingOutputMode() == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
						&& GvcfMathUtils.phredScaleErrorRate(PofF) >= options.STANDARD_CONFIDENCE_FOR_CALLING);
	}

	protected final boolean passesCallThreshold(double conf) {
		return conf >= options.STANDARD_CONFIDENCE_FOR_CALLING;
	}

	protected boolean forceKeepAllele(final Allele allele) {
		return options.getGenotypingOutputMode() == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
				|| options.ANNOTATE_ALL_SITES_WITH_PL;
	}

	protected final boolean passesEmitThreshold(double conf, boolean bestGuessIsRef) {
		return (options.getOutputMode() == OutputMode.EMIT_ALL_CONFIDENT_SITES || !bestGuessIsRef)
				&& conf >= options.STANDARD_CONFIDENCE_FOR_EMITTING;
	}

	private OutputAlleleSubset calculateOutputAlleleSubset(final AFCalculationResult afcr, final VariantContext vc) {
		final List<Allele> alleles = afcr.getAllelesUsedInGenotyping();

		final int alternativeAlleleCount = alleles.size() - 1;
		final Allele[] outputAlleles = new Allele[alternativeAlleleCount];
		final int[] mleCounts = new int[alternativeAlleleCount];
		int outputAlleleCount = 0;
		boolean siteIsMonomorphic = true;
		int referenceAlleleSize = 0;
		for (final Allele allele : alleles) {
			if (allele.isReference()) {
				referenceAlleleSize = allele.length();
			} else {
				// we want to keep the NON_REF symbolic allele but only in the
				// absence of a non-symbolic allele, e.g.
				// if we combined a ref / NON_REF gVCF with a ref / alt gVCF
				final boolean isNonRefWhichIsLoneAltAllele = alternativeAlleleCount == 1
						&& allele.equals(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE);
				final boolean isPlausible = afcr.isPolymorphicPhredScaledQual(allele,
						options.STANDARD_CONFIDENCE_FOR_CALLING);

				siteIsMonomorphic &= !isPlausible;
				boolean toOutput = (isPlausible || forceKeepAllele(allele) || isNonRefWhichIsLoneAltAllele);
				if (allele.equals(GaeaVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED)
						|| allele.equals(Allele.SPAN_DEL)) {
					toOutput &= coveredByDeletion(vc);
				}
				if (toOutput) {
					outputAlleles[outputAlleleCount] = allele;
					mleCounts[outputAlleleCount++] = afcr.getAlleleCountAtMLE(allele);
					recordDeletion(referenceAlleleSize, allele, vc);
				}
			}
		}

		return new OutputAlleleSubset(outputAlleleCount, outputAlleles, mleCounts, siteIsMonomorphic);
	}

	private void recordDeletion(final int referenceAlleleSize, final Allele allele, final VariantContext vc) {
		final int deletionSize = referenceAlleleSize - allele.length();

		// Allele ia a deletion
		if (deletionSize > 0) {
			final GenomeLocation genomeLoc = genomeLocParser.createGenomeLocation(vc.getContig(), vc.getStart(),
					vc.getStart() + deletionSize);
			upstreamDeletionsLoc.add(genomeLoc);
		}
	}

	private boolean coveredByDeletion(final VariantContext vc) {
		for (Iterator<GenomeLocation> it = upstreamDeletionsLoc.iterator(); it.hasNext();) {
			final GenomeLocation loc = it.next();
			if (!loc.getContig().equals(vc.getContig())) { // past contig
															// deletion.
				it.remove();
			} else if (loc.getStop() < vc.getStart()) { // past position in
														// current contig
														// deletion.
				it.remove();
			} else if (loc.getStart() == vc.getStart()) {
				// ignore this deletion, the symbolic one does not make
				// reference to it.
			} else { // deletion covers.
				return true;
			}
		}

		return false;
	}

	protected VariantCallContext calculateGenotypes(final VariantContext vc, final Model model,
			final boolean inheritAttributesFromInputVC, final boolean doAlleleSpecificCalcs) {

		// if input VC can't be genotyped, exit with either null VCC or, in case
		// where we need to emit all sites, an empty call
		SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
		if (hasTooManyAlternativeAlleles(vc) || vc.getNSamples() == 0)
			return null;

		final int defaultPloidy = options.getSamplePloidy();
		final int maxAltAlleles = options.getMaxAlternateAllele();
		final int maxNumPLValues = options.getMaxNumberPLValues();
		final AFCalculator afCalculatorForAlleleSubsetting = afCalculatorProvider.getInstance(vc,defaultPloidy,maxAltAlleles).setMaxNumPLValues(maxNumPLValues);
        final AFCalculator afCalculatorForQualScore = options.USE_NEW_AF_CALCULATOR ? newAFCalculator : afCalculatorForAlleleSubsetting;
		//慢在这里
        final AFCalculationResult AFresult = afCalculatorForQualScore.getLog10PNonRef(vc, defaultPloidy,maxAltAlleles, getAlleleFrequencyPriors(vc,defaultPloidy,model));
        final OutputAlleleSubset outputAlternativeAlleles = calculateOutputAlleleSubset(AFresult, vc);
		final double PoFGT0 = Math.pow(10, AFresult.getLog10PosteriorOfAFGT0());
		// note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't
		// nice
		final double log10Confidence = !outputAlternativeAlleles.siteIsMonomorphic
				|| options.getGenotypingOutputMode() == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
				|| options.ANNOTATE_ALL_SITES_WITH_PL ? AFresult.getLog10PosteriorOfAFEq0() + 0.0
						: AFresult.getLog10PosteriorOfAFGT0() + 0.0;

		// Add 0.0 removes -0.0 occurrences.
		final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;
		// return a null call if we don't pass the confidence cutoff or the most
		// likely allele frequency is zero
		// skip this if we are already looking at a vc with a NON_REF allele
		// i.e. if we are in GenotypeGVCFs
		if (!passesEmitThreshold(phredScaledConfidence, outputAlternativeAlleles.siteIsMonomorphic)
				&& !forceSiteEmission()
				&& outputAlternativeAlleles.alleles[0] != GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE) {
			return null;
		}

		// start constructing the resulting VC
		final GenomeLocationParser genomeLocParser = this.genomeLocParser;
		if (genomeLocParser == null)
			throw new IllegalStateException(
					"this UG engine was created without a valid genomeLocParser and no refContext was provided");
		final GenomeLocation loc = genomeLocParser.createGenomeLocation(vc);
		final List<Allele> outputAlleles = outputAlternativeAlleles.outputAlleles(vc.getReference());
		final VariantContextBuilder builder = new VariantContextBuilder("UG_JOINT_CALL", loc.getContig(),
				loc.getStart(), loc.getStop(), outputAlleles);
		builder.log10PError(log10Confidence);
		if (!passesCallThreshold(phredScaledConfidence))
			builder.filter(GaeaVCFConstants.LOW_QUAL_FILTER_NAME);

		// create the genotypes

		//final GenotypesContext genotypes = afCalculator.subsetAlleles(vc, defaultPloidy, outputAlleles, true);
		final GenotypesContext genotypes = afCalculatorForAlleleSubsetting.subsetAlleles(vc, defaultPloidy, outputAlleles, true);
		builder.genotypes(genotypes);
		// *** note that calculating strand bias involves overwriting data
		// structures, so we do that last
		final Map<String, Object> attributes = composeCallAttributes(inheritAttributesFromInputVC, vc,
				outputAlternativeAlleles.alternativeAlleleMLECounts(), outputAlternativeAlleles.siteIsMonomorphic,
				AFresult, outputAlternativeAlleles.outputAlleles(vc.getReference()), genotypes, model,
				doAlleleSpecificCalcs);
		builder.attributes(attributes);

		VariantContext vcCall = builder.make();
		return new VariantCallContext(vcCall, confidentlyCalled(phredScaledConfidence, PoFGT0));
	}

	protected Map<String, Object> composeCallAttributes(final boolean inheritAttributesFromInputVC,
			final VariantContext vc, final List<Integer> alleleCountsofMLE, final boolean bestGuessIsRef,
			final AFCalculationResult AFresult, final List<Allele> allAllelesToUse, final GenotypesContext genotypes,
			final Model model, final boolean doAlleleSpecificCalcs) {
		final HashMap<String, Object> attributes = new HashMap<>();

		// inherit attributes from input vc if requested
		if (inheritAttributesFromInputVC)
			attributes.putAll(vc.getAttributes());

		// add the MLE AC and AF annotations
		if (alleleCountsofMLE.size() > 0) {
			attributes.put(GaeaVCFConstants.MLE_ALLELE_COUNT_KEY, alleleCountsofMLE);
			final ArrayList<Double> MLEfrequencies = calculateMLEAlleleFrequencies(alleleCountsofMLE, genotypes);
			attributes.put(GaeaVCFConstants.MLE_ALLELE_FREQUENCY_KEY, MLEfrequencies);
		}

		if (options.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED)
			attributes.put(GaeaVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, vc.getAlternateAlleles().size());

		// since we don't have access to the AnnotationEngine, this argument
		// tells us whether we need per-allele QUALs for downstream annotations
		if (doAlleleSpecificCalcs) {
			List<Double> perAlleleQuals = new ArrayList<>();
			// Per-allele quals are not calculated for biallelic sites
			if (AFresult.getAllelesUsedInGenotyping().size() > 2) {
				for (final Allele a : allAllelesToUse) {
					if (a.isNonReference())
						perAlleleQuals.add(AFresult.getLog10PosteriorOfAFEq0ForAllele(a));
				}
			} else {
				perAlleleQuals.add(AFresult.getLog10PosteriorOfAFEq0());
			}

			attributes.put(GaeaVCFConstants.AS_QUAL_KEY, perAlleleQuals);
		}

		return attributes;
	}

	private ArrayList<Double> calculateMLEAlleleFrequencies(List<Integer> alleleCountsofMLE,
			GenotypesContext genotypes) {
		int AN = 0;
		for (final Genotype g : genotypes)
			for (final Allele a : g.getAlleles())
				if (!a.isNoCall())
					AN++;

		final ArrayList<Double> MLEfrequencies = new ArrayList<>(alleleCountsofMLE.size());
		// the MLEAC is allowed to be larger than the AN (e.g. in the case of
		// all PLs being 0, the GT is ./. but the exact model may arbitrarily
		// choose an AC>1)
		for (final int AC : alleleCountsofMLE)
			MLEfrequencies.add(Math.min(1.0, (double) AC / (double) AN));
		return MLEfrequencies;
	}

	public Set<VCFInfoHeaderLine> getAppropriateVCFInfoHeaders() {
		final Set<VCFInfoHeaderLine> headerInfo = new HashSet<>();
		if (options.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED)
			headerInfo.add(GaeaVcfHeaderLines.getInfoLine(GaeaVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));
		return headerInfo;
	}
}
