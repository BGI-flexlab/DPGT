package org.bgi.flexlab.gaea.tools.haplotypecaller.math;

import org.apache.commons.math3.analysis.integration.gauss.BaseRuleFactory;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;

public class HermiteRuleFactory extends BaseRuleFactory<Double> {
    /** &pi;<sup>1/2</sup> */
    private static final double SQRT_PI = 1.77245385090551602729;
    /** &pi;<sup>-1/4</sup> */
    private static final double H0 = 7.5112554446494248286e-1;
    /** &pi;<sup>-1/4</sup> &radic;2 */
    private static final double H1 = 1.0622519320271969145;

    /** {@inheritDoc} */
    @Override
    protected Pair<Double[], Double[]> computeRule(int numberOfPoints)
        throws DimensionMismatchException {

        if (numberOfPoints == 1) {
            // Break recursion.
            return new Pair<Double[], Double[]>(new Double[] { 0d },
                                                new Double[] { SQRT_PI });
        }

        // Get previous rule.
        // If it has not been computed yet it will trigger a recursive call
        // to this method.
        final int lastNumPoints = numberOfPoints - 1;
        final Double[] previousPoints = getRuleInternal(lastNumPoints).getFirst();

        // Compute next rule.
        final Double[] points = new Double[numberOfPoints];
        final Double[] weights = new Double[numberOfPoints];

        final double sqrtTwoTimesLastNumPoints = FastMath.sqrt(2 * lastNumPoints);
        final double sqrtTwoTimesNumPoints = FastMath.sqrt(2 * numberOfPoints);

        // Find i-th root of H[n+1] by bracketing.
        final int iMax = numberOfPoints / 2;
        for (int i = 0; i < iMax; i++) {
            // Lower-bound of the interval.
            double a = (i == 0) ? -sqrtTwoTimesLastNumPoints : previousPoints[i - 1].doubleValue();
            // Upper-bound of the interval.
            double b = (iMax == 1) ? -0.5 : previousPoints[i].doubleValue();

            // H[j-1](a)
            double hma = H0;
            // H[j](a)
            double ha = H1 * a;
            // H[j-1](b)
            double hmb = H0;
            // H[j](b)
            double hb = H1 * b;
            for (int j = 1; j < numberOfPoints; j++) {
                // Compute H[j+1](a) and H[j+1](b)
                final double jp1 = j + 1;
                final double s = FastMath.sqrt(2 / jp1);
                final double sm = FastMath.sqrt(j / jp1);
                final double hpa = s * a * ha - sm * hma;
                final double hpb = s * b * hb - sm * hmb;
                hma = ha;
                ha = hpa;
                hmb = hb;
                hb = hpb;
            }

            // Now ha = H[n+1](a), and hma = H[n](a) (same holds for b).
            // Middle of the interval.
            double c = 0.5 * (a + b);
            // P[j-1](c)
            double hmc = H0;
            // P[j](c)
            double hc = H1 * c;
            boolean done = false;
            while (!done) {
                done = b - a <= Math.ulp(c);
                hmc = H0;
                hc = H1 * c;
                for (int j = 1; j < numberOfPoints; j++) {
                    // Compute H[j+1](c)
                    final double jp1 = j + 1;
                    final double s = FastMath.sqrt(2 / jp1);
                    final double sm = FastMath.sqrt(j / jp1);
                    final double hpc = s * c * hc - sm * hmc;
                    hmc = hc;
                    hc = hpc;
                }
                // Now h = H[n+1](c) and hm = H[n](c).
                if (!done) {
                    if (ha * hc < 0) {
                        b = c;
                        hmb = hmc;
                        hb = hc;
                    } else {
                        a = c;
                        hma = hmc;
                        ha = hc;
                    }
                    c = 0.5 * (a + b);
                }
            }
            final double d = sqrtTwoTimesNumPoints * hmc;
            final double w = 2 / (d * d);

            points[i] = c;
            weights[i] = w;

            final int idx = lastNumPoints - i;
            points[idx] = -c;
            weights[idx] = w;
        }

        // If "numberOfPoints" is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (numberOfPoints % 2 != 0) {
            double hm = H0;
            for (int j = 1; j < numberOfPoints; j += 2) {
                final double jp1 = j + 1;
                hm = -FastMath.sqrt(j / jp1) * hm;
            }
            final double d = sqrtTwoTimesNumPoints * hm;
            final double w = 2 / (d * d);

            points[iMax] = 0d;
            weights[iMax] = w;
        }

        return new Pair<Double[], Double[]>(points, weights);
    }
}

