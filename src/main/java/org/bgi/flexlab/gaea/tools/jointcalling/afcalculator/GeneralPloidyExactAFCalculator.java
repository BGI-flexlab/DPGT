package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.genotyer.genotypecaller.ExactACcounts;
import org.bgi.flexlab.gaea.tools.genotyer.genotypecaller.ExactACset;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeAlleleCounts;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculator;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculators;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.MathUtils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public class GeneralPloidyExactAFCalculator extends ExactAFCalculator{
	protected GeneralPloidyExactAFCalculator() {
    }

    @Override
    protected GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse) {
        return subsetAlleles(vc, defaultPloidy, allelesToUse, false);
    }

    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy, final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        combineSinglePools(vc.getGenotypes(), defaultPloidy, vc.getNAlleles(), log10AlleleFrequencyPriors);
        return getResultFromFinalState(vc, log10AlleleFrequencyPriors, stateTracker);
    }

    /**
     * Simple wrapper class to hold values of combined pool likelihoods.
     * For fast hashing and fast retrieval, there's a hash map that shadows main list.
     *
     */
    static class CombinedPoolLikelihoods {
        private LinkedList<ExactACset> alleleCountSetList;
        private HashMap<ExactACcounts,ExactACset> conformationMap;
        private double maxLikelihood;


        public CombinedPoolLikelihoods() {
            // final int numElements = GenotypeLikelihoods.numLikelihoods();
            alleleCountSetList = new LinkedList<>();
            conformationMap = new HashMap<>();
            maxLikelihood = Double.NEGATIVE_INFINITY;
        }

        public void add(ExactACset set) {
            alleleCountSetList.add(set);
            conformationMap.put(set.getACcounts(), set);
            final double likelihood = set.getLog10Likelihoods()[0];

            if (likelihood > maxLikelihood )
                maxLikelihood = likelihood;

        }

        public boolean hasConformation(int[] ac) {
            return conformationMap.containsKey(new ExactACcounts(ac));

        }

        public double getLikelihoodOfConformation(int[] ac) {
            return conformationMap.get(new ExactACcounts(ac)).getLog10Likelihoods()[0];
        }

        public double getGLOfACZero() {
            return alleleCountSetList.get(0).getLog10Likelihoods()[0]; // AC 0 is always at beginning of list
        }

        public int getLength() {
            return alleleCountSetList.size();
        }
    }


    @Override
    protected void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums) {
        final int numOriginalAltAlleles = likelihoodSums.length;
        final GenotypesContext genotypes = vc.getGenotypes();
        for ( final Genotype genotype : genotypes.iterateInSampleNameOrder() ) {
            if (!genotype.hasPL())
                continue;
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (MathUtils.sum(gls) >= GaeaGvcfVariantContextUtils.SUM_GL_THRESH_NOCALL)
                continue;

            final int PLindexOfBestGL = MathUtils.maxElementIndex(gls);

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

    /**
     * Simple non-optimized version that combines GLs from several pools and produces global AF distribution.
     * @param GLs                              Inputs genotypes context with per-pool GLs
     * @param numAlleles                       Number of alternate alleles
     * @param log10AlleleFrequencyPriors       Frequency priors
     */
    protected void combineSinglePools(final GenotypesContext GLs,
                                      final int defaultPloidy,
                                      final int numAlleles,
                                      final double[] log10AlleleFrequencyPriors) {

        // Combine each pool incrementally - likelihoods will be renormalized at each step

        // first element: zero ploidy, e.g. trivial degenerate distribution
        final int numAltAlleles = numAlleles - 1;
        final int[] zeroCounts = new int[numAlleles];
        final ExactACset set = new ExactACset(1, new ExactACcounts(zeroCounts));
        set.getLog10Likelihoods()[0] = 0.0;
        final StateTracker stateTracker = getStateTracker(false,numAltAlleles);
        int combinedPloidy = 0;
        CombinedPoolLikelihoods combinedPoolLikelihoods = new CombinedPoolLikelihoods();
        combinedPoolLikelihoods.add(set);

        for (final Genotype genotype : GLs.iterateInSampleNameOrder()) {
            // recover gls and check if they qualify.
            if (!genotype.hasPL())
                continue;
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (MathUtils.sum(gls) >= GaeaGvcfVariantContextUtils.SUM_GL_THRESH_NOCALL)
                continue;
            stateTracker.reset();
            final int declaredPloidy = genotype.getPloidy();
            final int ploidy = declaredPloidy < 1 ? defaultPloidy : declaredPloidy;
            // they do qualify so we proceed.
            combinedPoolLikelihoods = fastCombineMultiallelicPool(combinedPoolLikelihoods, gls,
                    combinedPloidy, ploidy, numAlleles, log10AlleleFrequencyPriors, stateTracker);
            combinedPloidy = ploidy + combinedPloidy; // total number of chromosomes in combinedLikelihoods
        }
        if (combinedPloidy == 0)
            stateTracker.setLog10LikelihoodOfAFzero(0.0);
    }

    private CombinedPoolLikelihoods fastCombineMultiallelicPool(final CombinedPoolLikelihoods originalPool,
                                                               double[] newGL,
                                                               int originalPloidy,
                                                               int newGLPloidy,
                                                               int numAlleles,
                                                               final double[] log10AlleleFrequencyPriors,
                                                               final StateTracker stateTracker) {
        final LinkedList<ExactACset> ACqueue = new LinkedList<>();
        // mapping of ExactACset indexes to the objects
        final HashMap<ExactACcounts, ExactACset> indexesToACset = new HashMap<>();
        final CombinedPoolLikelihoods newPool = new CombinedPoolLikelihoods();

        // add AC=0 to the queue
        final int[] zeroCounts = new int[numAlleles];
        final int newPloidy = originalPloidy + newGLPloidy;
        zeroCounts[0] = newPloidy;

        ExactACset zeroSet = new ExactACset(1, new ExactACcounts(zeroCounts));

        ACqueue.add(zeroSet);
        indexesToACset.put(zeroSet.getACcounts(), zeroSet);

        // keep processing while we have AC conformations that need to be calculated
        while ( !ACqueue.isEmpty() ) {
            stateTracker.incNEvaluations();
            // compute log10Likelihoods
            final ExactACset ACset = ACqueue.remove();

            calculateACConformationAndUpdateQueue(ACset, newPool, originalPool, newGL, log10AlleleFrequencyPriors, originalPloidy, newGLPloidy, ACqueue, indexesToACset, stateTracker);

            // clean up memory
            indexesToACset.remove(ACset.getACcounts());
        }
        return newPool;
    }

    // todo - refactor, function almost identical except for log10LofK computation in GeneralPloidyGenotypeLikelihoods
    /**
     *
     * @param set                       ExactACset holding conformation to be computed
     * @param newPool                   New pool likelihood holder
     * @param originalPool              Original likelihood holder
     * @param newGL                     New pool GL vector to combine
     * @param log10AlleleFrequencyPriors Prior object
     * @param originalPloidy             Total ploidy of original combined pool
     * @param newGLPloidy                Ploidy of GL vector
     * @param ACqueue                    Queue of conformations to compute
     * @param indexesToACset             AC indices of objects in queue
     * @return                           max log likelihood
     */
    private double calculateACConformationAndUpdateQueue(final ExactACset set,
                                                         final CombinedPoolLikelihoods newPool,
                                                         final CombinedPoolLikelihoods originalPool,
                                                         final double[] newGL,
                                                         final double[] log10AlleleFrequencyPriors,
                                                         final int originalPloidy,
                                                         final int newGLPloidy,
                                                         final LinkedList<ExactACset> ACqueue,
                                                         final HashMap<ExactACcounts, ExactACset> indexesToACset,
                                                         final StateTracker stateTracker) {

        // compute likelihood in "set" of new set based on original likelihoods
        final int numAlleles = set.getACcounts().getCounts().length;
        final int newPloidy = set.getACsum();
        final double log10LofK = computeLofK(set, originalPool, newGL, log10AlleleFrequencyPriors, numAlleles, originalPloidy, newGLPloidy, stateTracker);


        // add to new pool
        if (!Double.isInfinite(log10LofK))
            newPool.add(set);

        if ( stateTracker.abort(log10LofK, set.getACcounts(), true, true) )
            return log10LofK;

        // iterate over higher frequencies if possible
        // by convention, ACcounts contained in set have full vector of possible pool ac counts including ref count.
        // so, if first element is zero, it automatically means we have no wiggle since we're in a corner of the conformation space
        final int ACwiggle = set.getACcounts().getCounts()[0];
        if ( ACwiggle == 0 ) // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return log10LofK;


        // add conformations for other cases
        for ( int allele = 1; allele < numAlleles; allele++ ) {
            final int[] ACcountsClone = set.getACcounts().getCounts().clone();
            ACcountsClone[allele]++;
            // is this a valid conformation?
            int altSum = (int)MathUtils.sum(ACcountsClone) - ACcountsClone[0];
            ACcountsClone[0] = newPloidy - altSum;
            if (ACcountsClone[0] < 0)
                continue;


            GeneralPloidyGenotypeLikelihoods.updateACset(ACcountsClone, ACqueue, indexesToACset);
        }


        return log10LofK;
    }


    /**
     * Compute likelihood of a particular AC conformation and update AFresult object
     * @param set                     Set of AC counts to compute
     * @param firstGLs                  Original pool likelihoods before combining
     * @param secondGL                  New GL vector with additional pool
     * @param log10AlleleFrequencyPriors     Allele frequency priors
     * @param numAlleles                Number of alleles (including ref)
     * @param ploidy1                   Ploidy of original pool (combined)
     * @param ploidy2                   Ploidy of new pool
     * @return                          log-likelihood of requested conformation
     */
    private double computeLofK(final ExactACset set,
                               final CombinedPoolLikelihoods firstGLs,
                               final double[] secondGL,
                               final double[] log10AlleleFrequencyPriors,
                               final int numAlleles, final int ploidy1, final int ploidy2, final StateTracker stateTracker) {

        final int newPloidy = ploidy1 + ploidy2;

        // sanity check
        int totalAltK = set.getACsum();
        if (newPloidy != totalAltK)
            throw new UserException("BUG: inconsistent sizes of set.getACsum and passed ploidy values");

        totalAltK -= set.getACcounts().getCounts()[0];
        // totalAltK has sum of alt alleles of conformation now


        // special case for k = 0 over all k
        if ( totalAltK == 0 ) {   // all-ref case
            final double log10Lof0 = firstGLs.getGLOfACZero() + secondGL[HOM_REF_INDEX];
            set.getLog10Likelihoods()[0] = log10Lof0;
            stateTracker.setLog10LikelihoodOfAFzero(log10Lof0);
            stateTracker.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);
            return log10Lof0;

        }   else {

            // initialize result with denominator
            // ExactACset holds by convention the conformation of all alleles, and the sum of all allele count is just the ploidy.
            // To compute n!/k1!k2!k3!... we need to compute first n!/(k2!k3!...) and then further divide by k1! where k1=ploidy-sum_k_i

            int[] currentCount = set.getACcounts().getCounts();
            double denom =  -GvcfMathUtils.log10MultinomialCoefficient(newPloidy, currentCount);

            // for current conformation, get all possible ways to break vector K into two components G1 and G2
            final GeneralPloidyGenotypeLikelihoods.SumIterator innerIterator = new GeneralPloidyGenotypeLikelihoods.SumIterator(numAlleles,ploidy2);
            set.getLog10Likelihoods()[0] = Double.NEGATIVE_INFINITY;
            while (innerIterator.hasNext()) {
                // check if breaking current conformation into g1 and g2 is feasible.
                final int[] acCount2 = innerIterator.getCurrentVector();
                final int[] acCount1 = GvcfMathUtils.vectorDiff(currentCount, acCount2);
                final int idx2 = innerIterator.getLinearIndex();
                // see if conformation is valid and if original pool had this conformation
                // for conformation to be valid, all elements of g2 have to be <= elements of current AC set
                if (isValidConformation(acCount1,ploidy1) && firstGLs.hasConformation(acCount1)) {
                    final double gl2 = secondGL[idx2];
                    if (!Double.isInfinite(gl2)) {
                        final double firstGL = firstGLs.getLikelihoodOfConformation(acCount1);
                        final double num1 = GvcfMathUtils.log10MultinomialCoefficient(ploidy1, acCount1);
                        final double num2 = GvcfMathUtils.log10MultinomialCoefficient(ploidy2, acCount2);
                        final double sum = firstGL + gl2 + num1 + num2;

                        set.getLog10Likelihoods()[0] = MathUtils.approximateLog10SumLog10(set.getLog10Likelihoods()[0], sum);
                    }
                }
                innerIterator.next();
            }

            set.getLog10Likelihoods()[0] += denom;
        }

        double log10LofK = set.getLog10Likelihoods()[0];

        // update the MLE if necessary
        final int altCounts[] = Arrays.copyOfRange(set.getACcounts().getCounts(),1, set.getACcounts().getCounts().length);
        // TODO -- GUILLERMO THIS CODE MAY PRODUCE POSITIVE LIKELIHOODS OR -INFINITY
        stateTracker.updateMLEifNeeded(Math.max(log10LofK, -Double.MAX_VALUE), altCounts);

        // apply the priors over each alternate allele
        for (final int ACcount : altCounts ) {
            if ( ACcount > 0 )
                log10LofK += log10AlleleFrequencyPriors[ACcount];
        }
        // TODO -- GUILLERMO THIS CODE MAY PRODUCE POSITIVE LIKELIHOODS OR -INFINITY
        stateTracker.updateMAPifNeeded(Math.max(log10LofK, -Double.MAX_VALUE), altCounts);

        return log10LofK;
    }

    /**
     * Small helper routine - is a particular AC conformation vector valid? ie are all elements non-negative and sum to ploidy?
     * @param set                            AC conformation vector
     * @param ploidy                         Ploidy of set
     * @return                               Valid conformation
     */
    private static boolean isValidConformation(final int[] set, final int ploidy) {
        int sum=0;
        for (final int ac: set) {
            if (ac < 0)
                return false;
            sum += ac;

        }

        return (sum == ploidy);
    }

    /**
     * From a given variant context, extract a given subset of alleles, and update genotype context accordingly,
     * including updating the PLs, ADs and SACs, and assign genotypes accordingly
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param defaultPloidy                     ploidy to assume in case that {@code vc} does not contain that information
     *                                          for a sample.
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes                   true: assign hard genotypes, false: leave as no-call
     * @return                                  GenotypesContext with new PLs, SACs and AD.
     */
    @Override
    public GenotypesContext subsetAlleles(final VariantContext vc, final int defaultPloidy,
                                          final List<Allele> allelesToUse,
                                          final boolean assignGenotypes) {

        final GenotypesContext result = GenotypesContext.create();

        // Subset genotypes for each sample
        for (final Genotype g : vc.getGenotypes()) // If it really needs to process order by sample name do so.
            result.add(subsetGenotypeAlleles(g, allelesToUse, vc, defaultPloidy, assignGenotypes));
        return GaeaGvcfVariantContextUtils.fixADFromSubsettedAlleles(result, vc, allelesToUse);
    }

    /**
     * From a given genotype, extract a given subset of alleles and update genotype PLs and SACs.
     * @param g                                 genotype to subset
     * @param allelesToUse                      alleles to subset
     * @param vc                                variant context with alleles and genotypes
     * @param defaultPloidy                     ploidy to assume in case that {@code vc} does not contain that information for a sample.
     * @param assignGenotypes                   true: assign hard genotypes, false: leave as no-call
     * @return                                  Genotypes with new PLs and SACs
     */
    private Genotype subsetGenotypeAlleles(final Genotype g, final List<Allele> allelesToUse, final VariantContext vc, final int defaultPloidy,
                                           boolean assignGenotypes) {
        final int ploidy = g.getPloidy() <= 0 ? defaultPloidy : g.getPloidy();
        if (!g.hasLikelihoods())
            return GenotypeBuilder.create(g.getSampleName(),GaeaGvcfVariantContextUtils.noCallAlleles(ploidy));
        else {
            // subset likelihood alleles
            final double[] newLikelihoods = subsetLikelihoodAlleles(g, allelesToUse, vc, ploidy);
            if (MathUtils.sum(newLikelihoods) > GaeaGvcfVariantContextUtils.SUM_GL_THRESH_NOCALL)
                return GenotypeBuilder.create(g.getSampleName(), GaeaGvcfVariantContextUtils.noCallAlleles(ploidy));
            else  // just now we would care about newSACs
                return subsetGenotypeAllelesWithLikelihoods(g, allelesToUse, vc, ploidy, assignGenotypes, newLikelihoods);
        }
    }

    /**
     * From a given genotype, extract a given subset of alleles and return the new PLs
     * @param g                                 genotype to subset
     * @param allelesToUse                      alleles to subset
     * @param vc                                variant context with alleles and genotypes
     * @param ploidy                            number of chromosomes
     * @return                                  the subsetted PLs
     */
    private double[] subsetLikelihoodAlleles(final Genotype g, final List<Allele> allelesToUse, final VariantContext vc, final int ploidy){

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;

        // create the new likelihoods array from the alleles we are allowed to use
        final double[] originalLikelihoods = g.getLikelihoods().getAsVector();

        if ( numOriginalAltAlleles != numNewAltAlleles ) {
            // might need to re-normalize the new likelihoods
            return MathUtils.normalizeFromLog10(GeneralPloidyGenotypeLikelihoods.subsetToAlleles(originalLikelihoods, ploidy, vc.getAlleles(), allelesToUse),
                    false, true);
        }
        else
            return originalLikelihoods;
    }

    /**
     * From a given genotype, subset the PLs and SACs
     * @param g                                 genotype to subset
     * @param allelesToUse                      alleles to subset
     * @param vc                                variant context with alleles and genotypes
     * @param ploidy                            number of chromosomes
     * @param assignGenotypes                   true: assign hard genotypes, false: leave as no-call
     * @param newLikelihoods                    the PL values
     * @return genotype with the subsetted PLsL and SACs
     */
    private Genotype subsetGenotypeAllelesWithLikelihoods(final Genotype g, final List<Allele> allelesToUse, final VariantContext vc, int ploidy,
                                                          final boolean assignGenotypes, final double[] newLikelihoods) {

        final GenotypeBuilder gb = new GenotypeBuilder(g);
        final String sampleName = g.getSampleName();

        // add likelihoods
        gb.PL(newLikelihoods);

        // get and add subsetted SACs
        final int[] newSACs = subsetSACAlleles(g, allelesToUse, vc);
        if (newSACs != null)
            gb.attribute(GaeaVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, newSACs);
        if (assignGenotypes)
            assignGenotype(gb, vc, sampleName, newLikelihoods, allelesToUse, ploidy);
        else
            gb.alleles(GaeaGvcfVariantContextUtils.noCallAlleles(ploidy));

        return gb.make();
    }

    /**
     * From a given genotype, extract a given subset of alleles and return the new SACs
     * @param g                             genotype to subset
     * @param allelesToUse                  alleles to subset
     * @param vc                            variant context with alleles and genotypes
     * @return                              the subsetted SACs
     */
    private int[] subsetSACAlleles(final Genotype g, final List<Allele> allelesToUse, final VariantContext vc){

        if ( !g.hasExtendedAttribute(GaeaVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY) )
            return null;

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;
        final List<Integer> sacIndexesToUse = numOriginalAltAlleles == numNewAltAlleles ? null : GaeaGvcfVariantContextUtils.determineSACIndexesToUse(vc, allelesToUse);

        return GaeaGvcfVariantContextUtils.makeNewSACs(g, sacIndexesToUse);
    }

    /**
     * Assign genotypes (GTs) to the samples in the VariantContext greedily based on the PLs
     *
     * @param gb                   the GenotypeBuilder to modify
     * @param vc                   the VariantContext
     * @param sampleName           the sample name
     * @param newLikelihoods       the PL array
     * @param allelesToUse         the list of alleles to choose from (corresponding to the PLs)
     * @param numChromosomes       Number of chromosomes per pool
     */
    private void assignGenotype(final GenotypeBuilder gb,
                                final VariantContext vc,
                                final String sampleName,
                                final double[] newLikelihoods,
                                final List<Allele> allelesToUse,
                                final int numChromosomes) {
        final int numNewAltAlleles = allelesToUse.size() - 1;

        // find the genotype with maximum likelihoods
        final int PLindex = numNewAltAlleles == 0 ? 0 : MathUtils.maxElementIndex(newLikelihoods);
        final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(numChromosomes,allelesToUse.size());
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(PLindex);

        gb.alleles(alleleCounts.asAlleleList(allelesToUse));

        removePLsIfMaxNumPLValuesExceeded(gb, vc, sampleName, newLikelihoods);

        // TODO - deprecated so what is the appropriate method to call?
        if ( numNewAltAlleles > 0 )
            gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
    }
}
