package org.bgi.flexlab.gaea.tools.haplotypecaller;

import java.util.stream.DoubleStream;

import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.ParamUtils;
import org.bgi.flexlab.gaea.util.Utils;

public final class RefVsAnyResult {
    /**
     * The genotype likelihoods for ref/ref ref/non-ref non-ref/non-ref
     */
    private final double[] genotypeLikelihoods;

    private int refDepth = 0;
    private int nonRefDepth = 0;

    /**
     * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
     * @param likelihoodCapacity the required capacity of the likelihood array, should match the possible number of
     *                           genotypes given the number of alleles (always 2), ploidy (arbitrary) less the genotyping
     *                           model non-sense genotype count if applies.
     * @throws IllegalArgumentException if {@code likelihoodCapacity} is negative.
     */
    public RefVsAnyResult(final int likelihoodCapacity) {
        ParamUtils.isPositiveOrZero(likelihoodCapacity, "likelihood capacity is negative");
        genotypeLikelihoods = new double[likelihoodCapacity];
    }

    public double[] getGenotypeLikelihoods() {
        return genotypeLikelihoods;
    }

    /**
     * @return Get the DP (sum of AD values)
     */
    public int getDP() { return refDepth + nonRefDepth; }

    /**
     * Return the AD fields. Returns a newly allocated array every time.
     */
    public int[] getAD(){
        return new int[]{refDepth, nonRefDepth};
    }

    /**
     * Returns (a copy of) the array of genotype likelihoods
     * Caps the het and hom var likelihood values by the hom ref likelihood.
     * The capping is done on the fly.
     */
    public double[] getGenotypeLikelihoodsCappedByHomRefLikelihood(){
        return DoubleStream.of(genotypeLikelihoods).map(d -> Math.min(d, genotypeLikelihoods[0])).toArray();
    }

    public void incrementRefAD(final int by){
        Utils.validateArg(by >= 0, "expected a non-negative number but got " + by);
        refDepth += by;
    }

    public void incrementNonRefAD(final int by){
        Utils.validateArg(by >= 0, "expected a non-negative number but got " + by);
        nonRefDepth += by;
    }

    public void addGenotypeLikelihood(final int idx, final double by) {
        genotypeLikelihoods[idx] += by;
    }
}
