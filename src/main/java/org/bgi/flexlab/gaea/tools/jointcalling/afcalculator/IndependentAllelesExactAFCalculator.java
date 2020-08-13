package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeAlleleCounts;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculator;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculators;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class IndependentAllelesExactAFCalculator extends ExactAFCalculator {
	/**
     * Array that caches the allele list that corresponds to the ith ploidy.
     *
     * <p>
     * Each position of the array <i>i</i >makes reference to a list that contains <i>i</i> copies of {@link Allele#NO_CALL}.
     * </p>
     *
     *<p>
     *     This array must be queried using {@link #biallelicNoCall(int)}, which will extend the cache
     *     to larger ploidies if needed.
     * </p>
     */
    private static volatile List<Allele>[] BIALLELIC_NOCALL = initialBiallelicNoCall(10);

    /**
     * Array that caches the allele list that corresponds to the ith ploidy.
     *
     * <p>
     * Each position of the array <i>i</i >makes reference to an array that contains
     * all-zero likelihoods with the number of genotypes that correspond
     * to a biallelic variant with ploidy <i>i</i>.
     * </p>
     *
     * <p>
     *     This array must be queried using {@link #biallelicNonInformativePls(int)}, which will extend the cache
     *     to larger ploidies if needed.
     * </p>
     */
    private static volatile int[][] BIALLELIC_NON_INFORMATIVE_PLS_BY_PLOIDY = initialBiallelicNonInformativePLsByPloidy(10);

    private static final Comparator<AFCalculationResult> AFCALC_RESULT_BY_PNONREF_COMPARATOR = new Comparator<AFCalculationResult>() {
        @Override
        public int compare(final AFCalculationResult o1, final AFCalculationResult o2) {
            return -1 * Double.compare(o1.getLog10PosteriorOfAFGT0(), o2.getLog10PosteriorOfAFGT0());
        }
    };

    private final ExactAFCalculator biallelicExactAFCalculator;

    protected IndependentAllelesExactAFCalculator(final ExactAFCalculator biallelicExactAFCalculator) {
        if (biallelicExactAFCalculator == null)
            throw new IllegalArgumentException("the biallelic exact AF calculator cannot be null");
        this.biallelicExactAFCalculator = biallelicExactAFCalculator;
    }

    /**
     * Creates a new calculator that delegates on {@link GeneralPloidyExactAFCalculator} to run
     * the exact model per allele.
     *
     * <p>
     *    Note: this constructor may be called using reflexion.
     * </p>
     */
    @SuppressWarnings("unused")
    protected IndependentAllelesExactAFCalculator() {
        this(new GeneralPloidyExactAFCalculator());
    }

    @Override
    protected void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums) {
        final int numOriginalAltAlleles = likelihoodSums.length;
        final GenotypesContext genotypes = vc.getGenotypes();
        for ( final Genotype genotype : genotypes.iterateInSampleNameOrder() ) {
            if (!genotype.hasPL())
                continue;
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (GvcfMathUtils.sum(gls) >= GaeaGvcfVariantContextUtils.SUM_GL_THRESH_NOCALL)
                continue;

            final int PLindexOfBestGL = GvcfMathUtils.maxElementIndex(gls);

            final double bestToHomRefDiffGL = PLindexOfBestGL == PL_INDEX_OF_HOM_REF ? 0.0 : gls[PLindexOfBestGL] - gls[PL_INDEX_OF_HOM_REF];
            final int declaredPloidy = genotype.getPloidy();
            final int ploidy = declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;

            final int[] acCount = GeneralPloidyGenotypeLikelihoods.getAlleleCountFromPLIndex(1 + numOriginalAltAlleles, ploidy, PLindexOfBestGL);
            // by convention, first count coming from getAlleleCountFromPLIndex comes from reference allele
            for (int k=1; k < acCount.length;k++)
                if (acCount[k] > 0 )
                    likelihoodSums[k-1].sum += acCount[k] * bestToHomRefDiffGL;
        }
    }

    @Override
    protected GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse) {
        return subsetAlleles(vc,defaultPloidy,allelesToUse,false);
    }

    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy, final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
    	    final List<AFCalculationResult> independentResultTrackers = computeAlleleIndependentExact(vc, defaultPloidy, 
    		log10AlleleFrequencyPriors);

        // Paranoia check:
        if ( independentResultTrackers.size() <= 1 )
            throw new IllegalStateException("Independent alleles model returned an empty list of results at VC " + vc);
        else if ( independentResultTrackers.size() == 2 ) {
            // fast path for the very common bi-allelic use case
            return independentResultTrackers.get(1);
        } else {
            final List<AFCalculationResult> alternativesOnly = new ArrayList<>(independentResultTrackers.size() - 1);
            for (int i = 1; i < independentResultTrackers.size(); i++)
                alternativesOnly.add(independentResultTrackers.get(i));
            // we are a multi-allelic, so we need to actually combine the results
            final List<AFCalculationResult> withMultiAllelicPriors = applyMultiAllelicPriors(alternativesOnly);
            return combineIndependentPNonRefs(vc, withMultiAllelicPriors, independentResultTrackers.get(0));
        }
    }

    protected final List<AFCalculationResult> applyMultiAllelicPriors(final List<AFCalculationResult> conditionalPNonRefResults) {
        final ArrayList<AFCalculationResult> sorted = new ArrayList<AFCalculationResult>(conditionalPNonRefResults);

        // sort the results, so the most likely allele is first
        Collections.sort(sorted, AFCALC_RESULT_BY_PNONREF_COMPARATOR);

        double lastPosteriorGt0 = sorted.get(0).getLog10PosteriorOfAFGT0();
        final double log10SingleAllelePriorOfAFGt0 = conditionalPNonRefResults.get(0).getLog10PriorOfAFGT0();

        for ( int i = 0; i < sorted.size(); i++ ) {
            if ( sorted.get(i).getLog10PosteriorOfAFGT0() > lastPosteriorGt0 )
                throw new IllegalStateException("pNonRefResults not sorted: lastPosteriorGt0 " + lastPosteriorGt0 + " but current is " + sorted.get(i).getLog10PosteriorOfAFGT0());

            final double log10PriorAFGt0 = (i + 1) * log10SingleAllelePriorOfAFGt0;
            final double log10PriorAFEq0 = Math.log10(1 - Math.pow(10, log10PriorAFGt0));
            final double[] thetaTONPriors = new double[] { log10PriorAFEq0, log10PriorAFGt0 };

            // bind pNonRef for allele to the posterior value of the AF > 0 with the new adjusted prior
            sorted.set(i, sorted.get(i).withNewPriors(GvcfMathUtils.normalizeFromLog10(thetaTONPriors, true)));
        }

        return sorted;
    }

    /**
     * Take the independent estimates of pNonRef for each alt allele and combine them into a single result
     *
     * Given n independent calculations for each of n alternate alleles create a single
     * combined AFCalcResult with:
     *
     * priors for AF == 0 equal to theta^N for the nth least likely allele
     * posteriors that reflect the combined chance that any alleles are segregating and corresponding
     * likelihoods
     * combined MLEs in the order of the alt alleles in vc
     *
     * @param sortedResultsWithThetaNPriors the pNonRef result for each allele independently
     */
    protected AFCalculationResult combineIndependentPNonRefs(final VariantContext vc,
                                                             final List<AFCalculationResult> sortedResultsWithThetaNPriors,
                                                             final AFCalculationResult combinedAltAllelesResult) {


        int nEvaluations = 0;
        final int nAltAlleles = sortedResultsWithThetaNPriors.size();
        final int[] alleleCountsOfMLE = new int[nAltAlleles];
        final Map<Allele, Double> log10pRefByAllele = new HashMap<>(nAltAlleles);

        // the sum of the log10 posteriors for AF == 0 and AF > 0 to determine joint probs

        for ( final AFCalculationResult sortedResultWithThetaNPriors : sortedResultsWithThetaNPriors ) {
            final Allele altAllele = sortedResultWithThetaNPriors.getAllelesUsedInGenotyping().get(1);
            final int altI = vc.getAlleles().indexOf(altAllele) - 1;

            // MLE of altI allele is simply the MLE of this allele in altAlleles
            alleleCountsOfMLE[altI] = sortedResultWithThetaNPriors.getAlleleCountAtMLE(altAllele);

            // bind pNonRef for allele to the posterior value of the AF > 0 with the new adjusted prior
            log10pRefByAllele.put(altAllele, sortedResultWithThetaNPriors.getLog10PosteriorOfAFEq0());

            // trivial -- update the number of evaluations
            nEvaluations += sortedResultWithThetaNPriors.nEvaluations;
        }

        return new IndependentAlleleAFCalculationResult(alleleCountsOfMLE, nEvaluations, vc.getAlleles(),
                // necessary to ensure all values < 0
        		GvcfMathUtils.normalizeFromLog10(new double[] { combinedAltAllelesResult.getLog10LikelihoodOfAFEq0(), combinedAltAllelesResult.getLog10LikelihoodOfAFGT0() }, true),
                // priors incorporate multiple alt alleles, must be normalized
        		GvcfMathUtils.normalizeFromLog10(new double[] { combinedAltAllelesResult.getLog10PriorOfAFEq0(), combinedAltAllelesResult.getLog10PriorOfAFGT0() }, true),
                log10pRefByAllele, sortedResultsWithThetaNPriors);
    }

    /**
     * Compute the conditional exact AFCalcResult for each allele in vc independently, returning
     * the result of each, in order of the alt alleles in VC
     *
     * @param vc the VariantContext we want to analyze, with at least 1 alt allele
     * @param log10AlleleFrequencyPriors the priors
     * @return a list of the AFCalcResults for each bi-allelic sub context of vc
     */
    protected final List<AFCalculationResult> computeAlleleIndependentExact(final VariantContext vc, final int defaultPloidy,
                                                                            final double[] log10AlleleFrequencyPriors) {
        final List<AFCalculationResult> results = new LinkedList<>();

        for ( final VariantContext subvc : makeAlleleConditionalContexts(vc, defaultPloidy) ) {
            final AFCalculationResult resultTracker = biallelicExactAFCalculator.getLog10PNonRef(subvc, defaultPloidy, vc.getNAlleles() - 1, log10AlleleFrequencyPriors);
            results.add(resultTracker);
        }

        return results;
    }

    /**
     * Returns the bi-allelic variant context for each alt allele in vc with bi-allelic likelihoods, in order
     *
     * @param vc the variant context to split.  Must have n.alt.alleles > 1
     * @return a bi-allelic variant context for each alt allele in vc
     */
    protected final List<VariantContext> makeAlleleConditionalContexts(final VariantContext vc, final int defaultPloidy) {
        final int nAlleles = vc.getNAlleles();

        // go through the work of ripping up the VC into its biallelic components
        final List<VariantContext> vcs = new LinkedList<>();

        for ( int alleleIndex = 0; alleleIndex < nAlleles; alleleIndex++ ) {
          vcs.add(biallelicCombinedGLs(vc, defaultPloidy, alleleIndex));
        }
        return vcs;
    }

    /**
     * Create a single bi-allelic variant context from rootVC with alt allele with index altAlleleIndex
     *
     * @param rootVC the root (potentially multi-allelic) variant context
     * @param alleleIndex index of the alt allele, from 0 == reference
     * @return a bi-allelic variant context based on rootVC
     */
    protected final VariantContext biallelicCombinedGLs(final VariantContext rootVC, final int defaultPloidy, final int alleleIndex) {
        if ( rootVC.isBiallelic() ) {
            return rootVC;
        } else {
            final int nAlleles = rootVC.getNAlleles();
            final List<Genotype> biallelicGenotypes = new ArrayList<>(rootVC.getNSamples());
            for ( final Genotype g : rootVC.getGenotypes() )
                biallelicGenotypes.add(combineGLs(g, defaultPloidy, alleleIndex, nAlleles));

            final VariantContextBuilder vcb = new VariantContextBuilder(rootVC);
            final Allele allele = alleleIndex == 0 ? rootVC.getReference() : rootVC.getAlternateAllele(alleleIndex - 1);
            vcb.alleles(alleleIndex == 0  ? Arrays.asList(allele, GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE) : Arrays.asList(rootVC.getReference(), allele));
            vcb.genotypes(biallelicGenotypes);
            return vcb.make();
        }
    }

    /**
     * Returns a new Genotype with the PLs of the multi-allelic original reduced to a bi-allelic case.
     *
     * <p>Uses the log-sum-exp trick in order to work well with very low PLs</p>
     *
     * <p>This is handled in the following way:</p>
     *
     * <p>Suppose we have for a A/B/C site the following GLs:</p>
     *
     * <p>AA AB BB AC BC CC</p>
     *
     * <p>and we want to get the bi-allelic GLs for X/B, where X is everything not B</p>
     *
     * <p>XX = AA + AC + CC (since X = A or C)<br/>
     * XB = AB + BC                           <br/>
     * BB = BB                                <br/>
     * </p>
     * <p>
     *     This implementation uses the log-sum-exp trick in order to avoid numeric instability (underflow).
     * </p>
     *
     * @param original the original multi-allelic genotype
     * @param alleleIndex the index of the alt allele we wish to keep in the bialleic case -- with ref == 0
     * @param numberOfAlleles the total number of alleles (alternatives + the reference).
     * @return a new biallelic genotype with appropriate PLs
     */
    private Genotype combineGLs(final Genotype original, final int defaultPloidy, final int alleleIndex, final int numberOfAlleles ) {

        final int declaredPloidy = original.getPloidy();
        final int ploidy = declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;
        if ( original.isNonInformative() )
            return new GenotypeBuilder(original).PL(biallelicNonInformativePls(ploidy)).alleles(biallelicNoCall(ploidy)).make();

        final int[] pls = original.getPL();

        final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(ploidy, numberOfAlleles);
        final double[] newPLs = new double[ploidy + 1];
        Arrays.fill(newPLs, Double.NEGATIVE_INFINITY);
        for (int i = 0; i < pls.length; i++) {
            final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(i);
            final int alleleCount = alleleCounts.alleleCountFor(alleleIndex);
            final int newPLIndex = alleleIndex == 0 ? ploidy - alleleCount : alleleCount;
            newPLs[newPLIndex] = GvcfMathUtils.approximateLog10SumLog10(newPLs[newPLIndex], -.1 * pls[i]);
        }

        return new GenotypeBuilder(original).PL(newPLs).alleles(biallelicNoCall(ploidy)).make();
    }

    private static List<Allele>[] initialBiallelicNoCall(final int initialCapacity) {
        final List<Allele>[] result = new List[initialCapacity + 1];
        for (int i = 0; i < result.length; i++) {
            result[i] = GaeaGvcfVariantContextUtils.noCallAlleles(i);
        }
        return result;
    }

    private static int[][] initialBiallelicNonInformativePLsByPloidy(final int initialCapacity) {
        final int[][] result = new int[initialCapacity + 1][];
        for (int i = 0; i < result.length; i++)
            result[i] = new int[i]; // { 0, 0, 0 ... 0} is the actual uninformative PL array.
        return result;
    }

    /**
     * Returns a cached array of non-informative PLs (all 0) for a given ploidy.
     * <p>
     *     Calling code must never change its elements.
     * </p>
     * @param ploidy the required ploidy.
     * @return never {@code null}.
     */
    private static int[] biallelicNonInformativePls (final int ploidy) {
        if (ploidy >= BIALLELIC_NON_INFORMATIVE_PLS_BY_PLOIDY.length) {
            return enlargeIfNecessaryBiallelicNonInformativePlsByPloidyAndGet(ploidy);
        } else {
            return BIALLELIC_NON_INFORMATIVE_PLS_BY_PLOIDY[ploidy];
        }
    }

    /**
     * Thread-safe expansion of {@link #BIALLELIC_NON_INFORMATIVE_PLS_BY_PLOIDY}.
     * @param ploidy the requested ploidy.
     * @return the uninformative likelihoods array for the requested ploidy.
     */
    private static synchronized int[] enlargeIfNecessaryBiallelicNonInformativePlsByPloidyAndGet(final int ploidy) {
        if (ploidy >= BIALLELIC_NON_INFORMATIVE_PLS_BY_PLOIDY.length) {
            final int[][] newValue = Arrays.copyOf(BIALLELIC_NON_INFORMATIVE_PLS_BY_PLOIDY, ploidy * 2);
            for (int i = newValue.length - 1; i >= BIALLELIC_NON_INFORMATIVE_PLS_BY_PLOIDY.length; i--)
                newValue[i] = new int[i]; // { 0, 0, 0.. } is the actual uninformative PL array.
            BIALLELIC_NON_INFORMATIVE_PLS_BY_PLOIDY = newValue;
        }
        return BIALLELIC_NON_INFORMATIVE_PLS_BY_PLOIDY[ploidy];
    }

    /**
     * Returns a cached list of no-call alleles {@link Allele#NO_CALL} that correspond to a given ploidy.
     * <p>
     *     Calling code must never change its elements.
     * </p>
     * @param ploidy the required ploidy.
     * @return never {@code null}.
     */
    private static List<Allele> biallelicNoCall (final int ploidy) {
        if (ploidy >= BIALLELIC_NOCALL.length) {
            return enlargeIfNecessaryBiallelicNoCallAndGet(ploidy);
        } else {
            return BIALLELIC_NOCALL[ploidy];
        }
    }

    /**
     * Thread-safe expansion of {@link #BIALLELIC_NOCALL}.
     * @param ploidy the requested ploidy.
     * @return the no-call allele list for the requested ploidy.
     */
    private static synchronized List<Allele> enlargeIfNecessaryBiallelicNoCallAndGet(final int ploidy) {
        if (ploidy >= BIALLELIC_NOCALL.length) {
            final List<Allele>[] newValue = Arrays.copyOf(BIALLELIC_NOCALL, ploidy * 2);
            for (int i = newValue.length - 1; i >= BIALLELIC_NOCALL.length; i--)
                newValue[i] = GaeaGvcfVariantContextUtils.noCallAlleles(i);
            BIALLELIC_NOCALL = newValue;
        }
        return BIALLELIC_NOCALL[ploidy];
    }

    @Override
    public GenotypesContext subsetAlleles(VariantContext vc, int defaultPloidy, List<Allele> allelesToUse, boolean assignGenotypes) {
    	Utils.nonNull(vc, "vc cann't be null!");
    	Utils.nonNull(allelesToUse, "alleleToUse cann't be null!");
        // the genotypes with PLs
        final GenotypesContext oldGTs = vc.getGenotypes();

        // samples
        final List<String> sampleIndices = oldGTs.getSampleNamesOrderedByName();

        // the new genotypes to create
        final GenotypesContext newGTs = GenotypesContext.create();

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;


        // create the new genotypes
        for ( int k = 0; k < oldGTs.size(); k++ ) {
            final Genotype g = oldGTs.get(sampleIndices.get(k));
            final int declaredPloidy = g.getPloidy();
            final int ploidy = declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;
            if ( !g.hasLikelihoods() ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(),GaeaGvcfVariantContextUtils.noCallAlleles(ploidy)));
                continue;
            }

            // create the new likelihoods array from the alleles we are allowed to use
            final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
            double[] newLikelihoods;

            // Optimization: if # of new alt alleles = 0 (pure ref call), keep original likelihoods so we skip normalization
            // and subsetting
            if ( numOriginalAltAlleles == numNewAltAlleles || numNewAltAlleles == 0) {
                newLikelihoods = originalLikelihoods;
            } else {
                newLikelihoods = GeneralPloidyGenotypeLikelihoods.subsetToAlleles(originalLikelihoods, ploidy, vc.getAlleles(), allelesToUse);

                // might need to re-normalize
                newLikelihoods = GvcfMathUtils.normalizeFromLog10(newLikelihoods, false, true);
            }

            // if there is no mass on the (new) likelihoods, then just no-call the sample
            if ( GvcfMathUtils.sum(newLikelihoods) > GaeaGvcfVariantContextUtils.SUM_GL_THRESH_NOCALL ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), GaeaGvcfVariantContextUtils.noCallAlleles(ploidy)));
            } else {
                final GenotypeBuilder gb = new GenotypeBuilder(g);
                final String sampleName = g.getSampleName();

                if ( numNewAltAlleles == 0 )
                    gb.noPL();
                else
                    gb.PL(newLikelihoods);

                // if we weren't asked to assign a genotype, then just no-call the sample
                if ( !assignGenotypes || GvcfMathUtils.sum(newLikelihoods) > GaeaGvcfVariantContextUtils.SUM_GL_THRESH_NOCALL )
                    gb.alleles(GaeaGvcfVariantContextUtils.noCallAlleles(ploidy));
                else
                    assignGenotype(gb, vc, sampleName, newLikelihoods, allelesToUse, ploidy);
                newGTs.add(gb.make());
            }
        }

        return GaeaGvcfVariantContextUtils.fixADFromSubsettedAlleles(newGTs, vc, allelesToUse);
    }


    /**
     * Assign genotypes (GTs) to the samples in the VariantContext greedily based on the PLs
     *
     * @param gb                   the GenotypeBuilder to modify
     * @param vc                   the VariantContext
     * @param sampleName           the sample name
     * @param newLikelihoods       the PL array
     * @param allelesToUse         the list of alleles to choose from (corresponding to the PLs)
     * @param numChromosomes        Number of chromosomes per pool
     */
    private void assignGenotype(final GenotypeBuilder gb,
                                final VariantContext vc,
                                final String sampleName,
                                final double[] newLikelihoods,
                                final List<Allele> allelesToUse,
                                final int numChromosomes) {
        final int numNewAltAlleles = allelesToUse.size() - 1;

        // find the genotype with maximum likelihoods
        final int PLindex = numNewAltAlleles == 0 ? 0 : GvcfMathUtils.maxElementIndex(newLikelihoods);
        final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(numChromosomes,allelesToUse.size());
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(PLindex);

        gb.alleles(alleleCounts.asAlleleList(allelesToUse));

        removePLsIfMaxNumPLValuesExceeded(gb, vc, sampleName, newLikelihoods);

        if ( numNewAltAlleles > 0 )
            gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
    }
}
