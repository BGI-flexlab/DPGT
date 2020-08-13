package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculators;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class IndependentAllelesDiploidExactAFCalculator extends DiploidExactAFCalculator {
	private final static int[] BIALLELIC_NON_INFORMATIVE_PLS = new int[]{0,0,0};
    private final static List<Allele> BIALLELIC_NOCALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    /**
     * Sorts AFCalcResults by their posteriors of AF > 0, so the
     */
    private final static class CompareAFCalculatorResultsByPNonRef implements Comparator<AFCalculationResult> {
        @Override
        public int compare(AFCalculationResult o1, AFCalculationResult o2) {
            return -1 * Double.compare(o1.getLog10PosteriorOfAFGT0(), o2.getLog10PosteriorOfAFGT0());
        }
    }

    private final static CompareAFCalculatorResultsByPNonRef compareAFCalcResultsByPNonRef = new CompareAFCalculatorResultsByPNonRef();

    /**
     * The AFCalc model we are using to do the bi-allelic computation
     */
    final AFCalculator biAlleleExactModel;

    protected IndependentAllelesDiploidExactAFCalculator() {
        super();
        biAlleleExactModel = new ReferenceDiploidExactAFCalculator();
    }

    @Override
    public AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                            final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        final List<AFCalculationResult> independentResultTrackers = computeAlleleIndependentExact(vc, defaultPloidy, log10AlleleFrequencyPriors);

        if ( independentResultTrackers.size() == 0 )
            throw new IllegalStateException("Independent alleles model returned an empty list of results at VC " + vc);

        if ( independentResultTrackers.size() == 1 ) {
            // fast path for the very common bi-allelic use case
            return independentResultTrackers.get(0);
        } else {
            final AFCalculationResult combinedAltAllelesResult = combineAltAlleleIndependentExact(vc,defaultPloidy,log10AlleleFrequencyPriors);
            // we are a multi-allelic, so we need to actually combine the results
            final List<AFCalculationResult> withMultiAllelicPriors = applyMultiAllelicPriors(independentResultTrackers);
            return combineIndependentPNonRefs(vc, withMultiAllelicPriors, combinedAltAllelesResult);
        }
    }

    private AFCalculationResult combineAltAlleleIndependentExact(final VariantContext vc, int defaultPloidy, double[] log10AlleleFrequencyPriors) {
        final VariantContext combinedAltAllelesVariantContext = makeCombinedAltAllelesVariantContext(vc);
        final AFCalculationResult resultTracker = biAlleleExactModel.getLog10PNonRef(combinedAltAllelesVariantContext, defaultPloidy, vc.getNAlleles() - 1, log10AlleleFrequencyPriors);
        return resultTracker;
    }

    private VariantContext makeCombinedAltAllelesVariantContext(final VariantContext vc) {
        final int nAltAlleles = vc.getNAlleles() - 1;

        if ( nAltAlleles == 1 )
            return vc;
        else {
            final VariantContextBuilder vcb = new VariantContextBuilder(vc);
            final Allele reference = vcb.getAlleles().get(0);
            vcb.alleles(Arrays.asList(reference, GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE));
            final int genotypeCount = GenotypeLikelihoodCalculators.genotypeCount(2, vc.getNAlleles());
            final double[] hetLikelihoods = new double[vc.getNAlleles() - 1];
            final double[] homAltLikelihoods = new double[genotypeCount - hetLikelihoods.length - 1];
            final double[] newLikelihoods = new double[3];
            final List<Genotype> newGenotypes = new ArrayList<>(vc.getNSamples());
            for (final Genotype oldGenotype : vc.getGenotypes()) {
                final GenotypeBuilder gb = new GenotypeBuilder(oldGenotype);
                final List<Allele> oldAlleles = oldGenotype.getAlleles();
                if (oldAlleles != null) {
                    final List<Allele> newAlleles = new ArrayList<>(oldAlleles.size());
                    for (int i = 0; i < oldAlleles.size(); i++) {
                        final Allele oldAllele = oldAlleles.get(i);
                        if (oldAllele.isReference())
                            newAlleles.add(reference);
                        else if (oldAllele.isNoCall())
                            newAlleles.add(Allele.NO_CALL);
                        else
                            newAlleles.add(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE);
                    }
                    gb.alleles(newAlleles);
                }
                if (oldGenotype.isNonInformative())
                    gb.PL(BIALLELIC_NON_INFORMATIVE_PLS);
                else if (combineAltAlleleLikelihoods(oldGenotype, genotypeCount, newLikelihoods, hetLikelihoods, homAltLikelihoods))
                    gb.PL(newLikelihoods);
                newGenotypes.add(gb.make());
            }
            return vcb.genotypes(newGenotypes).make();
        }
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
    	Utils.nonNull(vc, "vc cann't be null!");
    	Utils.nonNull(log10AlleleFrequencyPriors, "array cann't be null!");
        final List<AFCalculationResult> results = new LinkedList<AFCalculationResult>();

        for ( final VariantContext subvc : makeAlleleConditionalContexts(vc) ) {
            final AFCalculationResult resultTracker = biAlleleExactModel.getLog10PNonRef(subvc, defaultPloidy, vc.getNAlleles() - 1, log10AlleleFrequencyPriors);
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
    protected final List<VariantContext> makeAlleleConditionalContexts(final VariantContext vc) {
    	Utils.nonNull(vc, "vc cann't be null!");
        final int nAltAlleles = vc.getNAlleles() - 1;

        if ( nAltAlleles == 1 ) {
            // fast path for bi-allelic case.
            return Collections.singletonList(vc);
        } else {
            // go through the work of ripping up the VC into its biallelic components
            final List<VariantContext> vcs = new LinkedList<VariantContext>();

            for ( int altI = 0; altI < nAltAlleles; altI++ ) {
                vcs.add(biallelicCombinedGLs(vc, altI + 1));
            }

            return vcs;
        }
    }

    /**
     * Create a single bi-allelic variant context from rootVC with alt allele with index altAlleleIndex
     *
     * @param rootVC the root (potentially multi-allelic) variant context
     * @param altAlleleIndex index of the alt allele, from 0 == first alt allele
     * @return a bi-allelic variant context based on rootVC
     */
    protected final VariantContext biallelicCombinedGLs(final VariantContext rootVC, final int altAlleleIndex) {
        if ( rootVC.isBiallelic() ) {
            return rootVC;
        } else {
            final int nAlts = rootVC.getNAlleles() - 1;
            final List<Genotype> biallelicGenotypes = new ArrayList<Genotype>(rootVC.getNSamples());
            for ( final Genotype g : rootVC.getGenotypes() )
                biallelicGenotypes.add(combineGLsPrecise(g, altAlleleIndex, nAlts));

            final VariantContextBuilder vcb = new VariantContextBuilder(rootVC);
            final Allele altAllele = rootVC.getAlternateAllele(altAlleleIndex - 1);
            vcb.alleles(Arrays.asList(rootVC.getReference(), altAllele));
            vcb.genotypes(biallelicGenotypes);
            return vcb.make();
        }
    }

    private static final double PHRED_2_LOG10_COEFF = -.1;

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
     * BB = BB                                     <br/>
     * </p>
     * <p>
     *     This implementation use the log sum trick in order to avoid numeric inestability.
     * </p>
     *
     * @param original the original multi-allelic genotype
     * @param altIndex the index of the alt allele we wish to keep in the bialleic case -- with ref == 0
     * @param nAlts the total number of alt alleles
     * @return a new biallelic genotype with appropriate PLs
     */
    protected Genotype combineGLsPrecise(final Genotype original, final int altIndex, final int nAlts ) {

        if ( original.isNonInformative() )
            return new GenotypeBuilder(original).PL(BIALLELIC_NON_INFORMATIVE_PLS).alleles(BIALLELIC_NOCALL).make();

        if ( altIndex < 1 || altIndex > nAlts ) throw new IllegalStateException("altIndex must be between 1 and nAlts " + nAlts);

        final int[] pls = original.getPL();

        final int nAlleles = nAlts + 1;

        final int plCount = pls.length;

        double BB = 0;
        final double[] XBvalues = new double[nAlleles - 1];
        final double[] XXvalues = new double[plCount - nAlleles];

        int xbOffset = 0;
        int xxOffset = 0;
        for ( int index = 0; index < plCount; index++ ) {
            final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair pair = GenotypeLikelihoods.getAllelePair(index);
            int i = pair.alleleIndex1;
            int j = pair.alleleIndex2;
            if (i == j) {
              if (i == altIndex) BB = PHRED_2_LOG10_COEFF * pls[index]; else XXvalues[xxOffset++] = PHRED_2_LOG10_COEFF * pls[index];
            } else if (i == altIndex || j == altIndex)
              XBvalues[xbOffset++] = PHRED_2_LOG10_COEFF * pls[index];
            else
              XXvalues[xxOffset++] = PHRED_2_LOG10_COEFF * pls[index];
        }

        final double XB = GvcfMathUtils.log10sumLog10(XBvalues);
        final double XX = GvcfMathUtils.log10sumLog10(XXvalues);

        final double[] GLs = new double[] { XX, XB, BB};
        return new GenotypeBuilder(original).PL(GLs).alleles(BIALLELIC_NOCALL).make();
    }

    protected final List<AFCalculationResult> applyMultiAllelicPriors(final List<AFCalculationResult> conditionalPNonRefResults) {
        final ArrayList<AFCalculationResult> sorted = new ArrayList<AFCalculationResult>(conditionalPNonRefResults);

        // sort the results, so the most likely allele is first
        Collections.sort(sorted, compareAFCalcResultsByPNonRef);

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
        final Map<Allele, Double> log10pRefByAllele = new HashMap<Allele, Double>(nAltAlleles);

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

    private boolean combineAltAlleleLikelihoods(final Genotype g, final int plMaxIndex, final double[] dest,
                                                final double[] hetLikelihoods, final double[] homAltLikelihoods) {

        final int[] pls = g.getPL();
        if (pls == null)
            return false;
        int hetNextIndex = 0;
        int homAltNextIndex = 0;
        for (int plIndex = 1; plIndex < plMaxIndex; plIndex++) {
            final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(plIndex);
            if (alleles.alleleIndex1 == 0 || alleles.alleleIndex2 == 0)
                hetLikelihoods[hetNextIndex++] = pls[plIndex] * PHRED_2_LOG10_COEFF;
            else
                homAltLikelihoods[homAltNextIndex++] = pls[plIndex] * PHRED_2_LOG10_COEFF;
        }
        dest[0] = pls[0] * PHRED_2_LOG10_COEFF;
        dest[1] = GvcfMathUtils.approximateLog10SumLog10(hetLikelihoods);
        dest[2] = GvcfMathUtils.approximateLog10SumLog10(homAltLikelihoods);
        return true;
    }
}
