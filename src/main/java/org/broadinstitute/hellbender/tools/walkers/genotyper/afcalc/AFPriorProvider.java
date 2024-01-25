package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import java.util.Arrays;

import org.broadinstitute.hellbender.utils.MathUtils;

/**
 * Class that produces allele-frequency priors.
 */
public abstract class AFPriorProvider {

    private double[][] priorByTotalPloidy;
    // log10 sum of prior[1:]
    private double[] log10sumByTotalploidy;

    protected AFPriorProvider() {

    }

    /**
     * Returns the priors given a total-ploidy (the total number of genome copies across all samples).
     *
     * <p>
     *     For performance sake the returned value is a direct reference ot the cached prior, so the client code
     *     must not modify its content.
     * </p>
     *
     * @param totalPloidy the requested total ploidy.
     *
     * @return never {@code null}, please do not modify its content. An array of {@code totalPloidy + 1} positions where
     *  the ith position is the log10(prior) probability of the an alternative allele AC to be exactly <i>i</i> copies in
     *  a draw of {@code totalPloidy} elements.
     */
    public double[] forTotalPloidy(final int totalPloidy) {
        if (totalPloidy < 0) {
            throw new IllegalArgumentException("the total-ploidy cannot be negative");
        }
        ensureCapacity(totalPloidy);
        final double[] cachedResult = priorByTotalPloidy[totalPloidy];
        if (cachedResult == null) {
            return priorByTotalPloidy[totalPloidy] = buildPriors(totalPloidy);
        } else {
            return cachedResult;
        }
    }


    /**
     * return log10(sum(priors[1:]))
     * For performance sake the returned value is a direct reference to the cached log10sum values.
     * @param totalPloidy the requested total ploidy.
     * @return log10(sum(priors[1:]))
     */
    public double getLog10sum(final int totalPloidy) {
        if (totalPloidy < 0) {
            throw new IllegalArgumentException("the total-ploidy cannot be negative");
        }
        ensureCapacity(totalPloidy);
        final double cachedResult = log10sumByTotalploidy[totalPloidy];
        if (cachedResult == 0.0) {
            final double[] cachedPrior = priorByTotalPloidy[totalPloidy];
            if (cachedPrior == null) {
                priorByTotalPloidy[totalPloidy] = buildPriors(totalPloidy);
            }
            log10sumByTotalploidy[totalPloidy] = MathUtils.log10sumLog10(priorByTotalPloidy[totalPloidy], 1);
            return log10sumByTotalploidy[totalPloidy];
        } else {
            return cachedResult;
        }
    }


    /**
     * Make sure that structures have enough capacity to hold the information up to the given total-ploidy.
     * @param totalPloidy
     */
    protected void ensureCapacity(final int totalPloidy) {
        if (totalPloidy < 0) {
            throw new IllegalArgumentException("the total-ploidy cannot be negative");
        }
        if (priorByTotalPloidy == null) {
            priorByTotalPloidy = new double[totalPloidy + 1][];  // just enough for those cases in where we have a fix total-ploidy.
            log10sumByTotalploidy = new double[totalPloidy + 1];
        } else if (priorByTotalPloidy.length - 1 < totalPloidy) {
            priorByTotalPloidy = Arrays.copyOf(priorByTotalPloidy, Math.max(priorByTotalPloidy.length << 1, totalPloidy + 1));
            log10sumByTotalploidy = Arrays.copyOf(log10sumByTotalploidy, Math.max(log10sumByTotalploidy.length << 1, totalPloidy + 1));
        }
    }

    /**
     * Given a total ploidy construct the allele prior probabilities array.
     * @param totalPloidy the target total-ploidy. Code can assume that is a non-negative number.
     *
     * @return never {@code null}, an array of exactly {@code totalPloidy + 1} positions that satisifed the
     *  contract {@link #forTotalPloidy(int) forTotalPloidy(totalPloidy)}.
     */
    protected abstract double[] buildPriors(final int totalPloidy);

}
