/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.util;

import org.apache.commons.math3.exception.NotFiniteNumberException;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.util.FastMath;
import org.bgi.flexlab.gaea.data.exception.UserException;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Collection;
import java.util.function.DoubleUnaryOperator;

public class MathUtils {

	public static DecimalFormat doubleformat = new DecimalFormat("0.000");
	public final static double LOG10_P_OF_ZERO = -1000000.0;
	public static final double LOG_ONE_HALF = -Math.log10(2.0);

	public static final double INV_SQRT_2_PI = 1.0 / Math.sqrt(2.0 * Math.PI);

	/**
	 * Private constructor. No instantiating this class!
	 */
	protected MathUtils() {
	}

	static {
		doubleformat.setRoundingMode(RoundingMode.HALF_UP);
	}

	public static final double[] log10Cache;
	public static final double[] log10FactorialCache;
	private static final int MAXN = 50000;
	private static final int LOG10_CACHE_SIZE = 4 * MAXN; // we need to be able
															// to go up to
															// 2*(2N) when
															// calculating some
															// of the
															// coefficients
	protected static final double JACOBIAN_LOG_TABLE_STEP = 0.001;
	protected static final double JACOBIAN_LOG_TABLE_INV_STEP = 1.0 / 0.001;
	protected static final double MAX_JACOBIAN_TOLERANCE = 8.0;
	protected static final int JACOBIAN_LOG_TABLE_SIZE = (int) (MAX_JACOBIAN_TOLERANCE / JACOBIAN_LOG_TABLE_STEP) + 1;
	protected static final double[] jacobianLogTable;
	
	public static final double FAIR_BINOMIAL_PROB_LOG10_0_5 = Math.log10(0.5);

	static {
		log10Cache = new double[LOG10_CACHE_SIZE];
		log10FactorialCache = new double[LOG10_CACHE_SIZE];
		jacobianLogTable = new double[JACOBIAN_LOG_TABLE_SIZE];

		log10Cache[0] = Double.NEGATIVE_INFINITY;
		for (int k = 1; k < LOG10_CACHE_SIZE; k++) {
			log10Cache[k] = Math.log10(k);
			log10FactorialCache[k] = log10FactorialCache[k - 1] + log10Cache[k];
		}

		for (int k = 0; k < JACOBIAN_LOG_TABLE_SIZE; k++) {
			jacobianLogTable[k] = Math.log10(1.0 + Math.pow(10.0, -((double) k) * JACOBIAN_LOG_TABLE_STEP));
		}
	}
	
	private static final class JacobianLogTable {
	    // if log(a) - log(b) > MAX_TOLERANCE, b is effectively treated as zero in approximateLogSumLog
	    // MAX_TOLERANCE = 8.0 introduces an error of at most one part in 10^8 in sums
	    public static final double MAX_TOLERANCE = 8.0;

	    //  Phred scores Q and Q+1 differ by 0.1 in their corresponding log-10 probabilities, and by
	    // 0.1 * log(10) in natural log probabilities.  Setting TABLE_STEP to an exact divisor of this
	    // quantity ensures that approximateSumLog in fact caches exact values for integer phred scores
	    private static final double TABLE_STEP = 0.0001;
	    private static final double INV_STEP = 1.0 / TABLE_STEP;
	    private static final double[] cache = new IndexRange(0, (int) (MAX_TOLERANCE / TABLE_STEP) + 1)
	            .mapToDouble(k -> Math.log10(1.0 + Math.pow(10.0, -k * TABLE_STEP)));

	    public static double get(final double difference) {
	        final int index = fastRound(difference * INV_STEP);
	        return cache[index];
	    }
	}

	public static boolean wellFormedDouble(final double val) {
		return !Double.isInfinite(val) && !Double.isNaN(val);
	}

	public static int fastRound(double d) {
		return (d > 0.0) ? (int) (d + 0.5d) : (int) (d - 0.5d);
	}

	public static int sum(byte[] x) {
		int total = 0;
		for (byte v : x)
			total += (int) v;
		return total;
	}

	public static double sum(double[] values) {
		double s = 0.0;
		for (double v : values)
			s += v;
		return s;
	}

	public static long sum(int[] x) {
		long total = 0;
		for (int v : x)
			total += v;
		return total;
	}

	public static int[] addArrays(int[] a, int[] b) {
		int[] c = new int[a.length];
		for (int i = 0; i < a.length; i++)
			c[i] = a[i] + b[i];
		return c;
	}

	public static byte compareDoubles(double a, double b) {
		return compareDoubles(a, b, 1e-6);
	}

	public static byte compareDoubles(double a, double b, double epsilon) {
		if (Math.abs(a - b) < epsilon) {
			return 0;
		}
		if (a > b) {
			return -1;
		}
		return 1;
	}

	public static double distanceSquared(final double[] x, final double[] y) {
		double dist = 0.0;
		for (int iii = 0; iii < x.length; iii++) {
			dist += (x[iii] - y[iii]) * (x[iii] - y[iii]);
		}
		return dist;
	}

	static {
		doubleformat.setRoundingMode(RoundingMode.HALF_UP);
	}

	public static int maxElementIndex(final double[] array) {
		return maxElementIndex(array, array.length);
	}

	public static int maxElementIndex(final double[] array, final int endIndex) {
		if (array == null || array.length == 0)
			throw new IllegalArgumentException("Array cannot be null!");

		int maxI = 0;
		for (int i = 1; i < endIndex; i++) {
			if (array[i] > array[maxI])
				maxI = i;
		}

		return maxI;
	}

	public static int minElementIndex(double[] array) {
		if (array == null || array.length == 0)
			throw new IllegalArgumentException("Array cannot be null!");

		int minI = 0;
		for (int i = 1; i < array.length; i++) {
			if (array[i] < array[minI])
				minI = i;
		}

		return minI;
	}

	public static double arrayMax(final double[] array) {
		return array[maxElementIndex(array)];
	}

	public static double arrayMax(final double[] array, final int endIndex) {
		return array[maxElementIndex(array, endIndex)];
	}

	public static double arrayMin(double[] array) {
		return array[minElementIndex(array)];
	}

	/**
	 * Converts LN to LOG10
	 *
	 * @param ln
	 *            log(x)
	 * @return log10(x)
	 */
	public static double lnToLog10(double ln) {
		return ln * Math.log10(Math.exp(1));
	}

	/**
	 * normalizes the log10-based array. ASSUMES THAT ALL ARRAY ENTRIES ARE <= 0
	 * (<= 1 IN REAL-SPACE).
	 *
	 * @param array
	 *            the array to be normalized
	 * @param takeLog10OfOutput
	 *            if true, the output will be transformed back into log10 units
	 * @return a newly allocated array corresponding the normalized values in
	 *         array, maybe log10 transformed
	 */
	public static double[] normalizeFromLog10(double[] array, boolean takeLog10OfOutput) {
		return normalizeFromLog10(array, takeLog10OfOutput, false);
	}

	public static double[] normalizeFromLog10(double[] array, boolean takeLog10OfOutput, boolean keepInLogSpace) {
		// for precision purposes, we need to add (or really subtract, since
		// they're
		// all negative) the largest value; also, we need to convert to
		// normal-space.
		double maxValue = arrayMax(array);

		// we may decide to just normalize in log space without converting to
		// linear space
		if (keepInLogSpace) {
			for (int i = 0; i < array.length; i++) {
				array[i] -= maxValue;
			}
			return array;
		}

		// default case: go to linear space
		double[] normalized = new double[array.length];

		for (int i = 0; i < array.length; i++)
			normalized[i] = Math.pow(10, array[i] - maxValue);

		// normalize
		double sum = 0.0;
		for (int i = 0; i < array.length; i++)
			sum += normalized[i];
		for (int i = 0; i < array.length; i++) {
			double x = normalized[i] / sum;
			if (takeLog10OfOutput) {
				x = Math.log10(x);
				if (x < LOG10_P_OF_ZERO || Double.isInfinite(x))
					x = array[i] - maxValue;
			}

			normalized[i] = x;
		}

		return normalized;
	}

	/**
	 * normalizes the log10-based array. ASSUMES THAT ALL ARRAY ENTRIES ARE <= 0
	 * (<= 1 IN REAL-SPACE).
	 *
	 * @param array
	 *            the array to be normalized
	 * @return a newly allocated array corresponding the normalized values in
	 *         array
	 */
	public static double[] normalizeFromLog10(double[] array) {
		return normalizeFromLog10(array, false);
	}

	public static double sumLog10(double[] log10values) {
		return Math.pow(10.0, log10sumLog10(log10values));
		// double s = 0.0;
		// for ( double v : log10values) s += Math.pow(10.0, v);
		// return s;
	}

	/**
	 * Same routine, unboxed types for efficiency
	 *
	 * @param x
	 *            First vector
	 * @param y
	 *            Second vector
	 * @return Vector of same length as x and y so that z[k] = x[k]+y[k]
	 */
	public static double[] vectorSum(double[] x, double[] y) {
		if (x.length != y.length)
			throw new UserException("BUG: Lengths of x and y must be the same");

		double[] result = new double[x.length];
		for (int k = 0; k < x.length; k++)
			result[k] = x[k] + y[k];

		return result;
	}

	public static double log10sumLog10(double[] log10values) {
		return log10sumLog10(log10values, 0);
	}

	public static double log10sumLog10(double[] log10p, int start) {
		return log10sumLog10(log10p, start, log10p.length);
	}

	public static double log10sumLog10(double[] log10p, int start, int finish) {
		double sum = 0.0;

		double maxValue = arrayMax(log10p, finish);
		if (maxValue == Double.NEGATIVE_INFINITY)
			return maxValue;

		for (int i = start; i < finish; i++) {
			sum += Math.pow(10.0, log10p[i] - maxValue);
		}

		return Math.log10(sum) + maxValue;
	}

	public static double approximateLog10SumLog10(double a, double b, double c) {
		return approximateLog10SumLog10(a, approximateLog10SumLog10(b, c));
	}

	public static double approximateLog10SumLog10(double small, double big) {
		// make sure small is really the smaller value
		if (small > big) {
			final double t = big;
			big = small;
			small = t;
		}

		if (small == Double.NEGATIVE_INFINITY || big == Double.NEGATIVE_INFINITY)
			return big;

		final double diff = big - small;
		if (diff >= MathUtils.MAX_JACOBIAN_TOLERANCE)
			return big;

		// OK, so |y-x| < tol: we use the following identity then:
		// we need to compute log10(10^x + 10^y)
		// By Jacobian logarithm identity, this is equal to
		// max(x,y) + log10(1+10^-abs(x-y))
		// we compute the second term as a table lookup with integer
		// quantization
		// we have pre-stored correction for 0,0.1,0.2,... 10.0
		final int ind = fastRound(diff * MathUtils.JACOBIAN_LOG_TABLE_INV_STEP); // hard
																					// rounding
		return big + MathUtils.jacobianLogTable[ind];
	}
	
	public static double approximateLog10SumLog10(final double[] vals) {
        return approximateLog10SumLog10(vals, vals.length);
    }
	
	public static double approximateLog10SumLog10(final double[] vals, final int endIndex) {
	    final int maxElementIndex = MathUtils.maxElementIndex(vals, endIndex);
	    double approxSum = vals[maxElementIndex];

	    for (int i = 0; i < endIndex; i++) {
	        if (i == maxElementIndex || vals[i] == Double.NEGATIVE_INFINITY) {
	            continue;
	        }

	        // if vals[i] isn't too tiny relative to the sum so far, add it; otherwise ignore it
	        final double diff = approxSum - vals[i];
	        approxSum += diff < JacobianLogTable.MAX_TOLERANCE ? JacobianLogTable.get(diff) : 0.0;
	    }

	    return approxSum;
	}
	
	public static int maxElementIndex(final double[] array, final int start, final int endIndex) {
	    Utils.nonNull(array, "array may not be null");
	    Utils.validateArg(array.length > 0, "array may not be empty");
	    Utils.validateArg(start <= endIndex, "Start cannot be after end.");

	    int maxI = start;
	    for (int i = (start+1); i < endIndex; i++) {
	        if (array[i] > array[maxI])
	            maxI = i;
	    }
	    return maxI;
	}
	
	public static double approximateLog10SumLog10(final double[] vals, final int fromIndex, final int toIndex) {
	    if (fromIndex == toIndex) return Double.NEGATIVE_INFINITY;
	    final int maxElementIndex = MathUtils.maxElementIndex(vals,fromIndex,toIndex);
	    double approxSum = vals[maxElementIndex];

	    for (int i = fromIndex; i < toIndex; i++) {
	        final double val;
	        if (i == maxElementIndex || (val = vals[i]) == Double.NEGATIVE_INFINITY)
	            continue;
	        final double diff = approxSum - val;
	        if (diff < JacobianLogTable.MAX_TOLERANCE)
	            approxSum += JacobianLogTable.get(diff);
	    }
	    return approxSum;
	}

	public static boolean goodLog10ProbVector(final double[] vector, final int expectedSize,
			final boolean shouldSumToOne) {
		if (vector.length != expectedSize)
			return false;

		for (final double pr : vector) {
			if (!goodLog10Probability(pr))
				return false;
		}

		if (shouldSumToOne && compareDoubles(sumLog10(vector), 1.0, 1e-4) != 0)
			return false;

		return true; // everything is good
	}

	public static boolean goodLog10Probability(final double result) {
		return result <= 0.0 && !Double.isInfinite(result) && !Double.isNaN(result);
	}

	/**
	 * Constants to simplify the log gamma function calculation.
	 */
	private static final double zero = 0.0, one = 1.0, half = .5, a0 = 7.72156649015328655494e-02,
			a1 = 3.22467033424113591611e-01, a2 = 6.73523010531292681824e-02, a3 = 2.05808084325167332806e-02,
			a4 = 7.38555086081402883957e-03, a5 = 2.89051383673415629091e-03, a6 = 1.19270763183362067845e-03,
			a7 = 5.10069792153511336608e-04, a8 = 2.20862790713908385557e-04, a9 = 1.08011567247583939954e-04,
			a10 = 2.52144565451257326939e-05, a11 = 4.48640949618915160150e-05, tc = 1.46163214496836224576e+00,
			tf = -1.21486290535849611461e-01, tt = -3.63867699703950536541e-18, t0 = 4.83836122723810047042e-01,
			t1 = -1.47587722994593911752e-01, t2 = 6.46249402391333854778e-02, t3 = -3.27885410759859649565e-02,
			t4 = 1.79706750811820387126e-02, t5 = -1.03142241298341437450e-02, t6 = 6.10053870246291332635e-03,
			t7 = -3.68452016781138256760e-03, t8 = 2.25964780900612472250e-03, t9 = -1.40346469989232843813e-03,
			t10 = 8.81081882437654011382e-04, t11 = -5.38595305356740546715e-04, t12 = 3.15632070903625950361e-04,
			t13 = -3.12754168375120860518e-04, t14 = 3.35529192635519073543e-04, u0 = -7.72156649015328655494e-02,
			u1 = 6.32827064025093366517e-01, u2 = 1.45492250137234768737e+00, u3 = 9.77717527963372745603e-01,
			u4 = 2.28963728064692451092e-01, u5 = 1.33810918536787660377e-02, v1 = 2.45597793713041134822e+00,
			v2 = 2.12848976379893395361e+00, v3 = 7.69285150456672783825e-01, v4 = 1.04222645593369134254e-01,
			v5 = 3.21709242282423911810e-03, s0 = -7.72156649015328655494e-02, s1 = 2.14982415960608852501e-01,
			s2 = 3.25778796408930981787e-01, s3 = 1.46350472652464452805e-01, s4 = 2.66422703033638609560e-02,
			s5 = 1.84028451407337715652e-03, s6 = 3.19475326584100867617e-05, r1 = 1.39200533467621045958e+00,
			r2 = 7.21935547567138069525e-01, r3 = 1.71933865632803078993e-01, r4 = 1.86459191715652901344e-02,
			r5 = 7.77942496381893596434e-04, r6 = 7.32668430744625636189e-06, w0 = 4.18938533204672725052e-01,
			w1 = 8.33333333333329678849e-02, w2 = -2.77777777728775536470e-03, w3 = 7.93650558643019558500e-04,
			w4 = -5.95187557450339963135e-04, w5 = 8.36339918996282139126e-04, w6 = -1.63092934096575273989e-03;

	/**
	 * Efficient rounding functions to simplify the log gamma function
	 * calculation double to long with 32 bit shift
	 */
	private static final int HI(double x) {
		return (int) (Double.doubleToLongBits(x) >> 32);
	}

	/**
	 * Efficient rounding functions to simplify the log gamma function
	 * calculation double to long without shift
	 */
	private static final int LO(double x) {
		return (int) Double.doubleToLongBits(x);
	}

	/**
	 * Most efficent implementation of the lnGamma (FDLIBM) Use via the
	 * log10Gamma wrapper method.
	 */
	private static double lnGamma(double x) {
		double t, y, z, p, p1, p2, p3, q, r, w;
		int i;

		int hx = HI(x);
		int lx = LO(x);

		/* purge off +-inf, NaN, +-0, and negative arguments */
		int ix = hx & 0x7fffffff;
		if (ix >= 0x7ff00000)
			return Double.POSITIVE_INFINITY;
		if ((ix | lx) == 0 || hx < 0)
			return Double.NaN;
		if (ix < 0x3b900000) { /* |x|<2**-70, return -log(|x|) */
			return -Math.log(x);
		}

		/* purge off 1 and 2 */
		if ((((ix - 0x3ff00000) | lx) == 0) || (((ix - 0x40000000) | lx) == 0))
			r = 0;
		/* for x < 2.0 */
		else if (ix < 0x40000000) {
			if (ix <= 0x3feccccc) { /* lgamma(x) = lgamma(x+1)-log(x) */
				r = -Math.log(x);
				if (ix >= 0x3FE76944) {
					y = one - x;
					i = 0;
				} else if (ix >= 0x3FCDA661) {
					y = x - (tc - one);
					i = 1;
				} else {
					y = x;
					i = 2;
				}
			} else {
				r = zero;
				if (ix >= 0x3FFBB4C3) {
					y = 2.0 - x;
					i = 0;
				} /* [1.7316,2] */
				else if (ix >= 0x3FF3B4C4) {
					y = x - tc;
					i = 1;
				} /* [1.23,1.73] */
				else {
					y = x - one;
					i = 2;
				}
			}

			switch (i) {
			case 0:
				z = y * y;
				p1 = a0 + z * (a2 + z * (a4 + z * (a6 + z * (a8 + z * a10))));
				p2 = z * (a1 + z * (a3 + z * (a5 + z * (a7 + z * (a9 + z * a11)))));
				p = y * p1 + p2;
				r += (p - 0.5 * y);
				break;
			case 1:
				z = y * y;
				w = z * y;
				p1 = t0 + w * (t3
						+ w * (t6 + w * (t9 + w * t12))); /* parallel comp */
				p2 = t1 + w * (t4 + w * (t7 + w * (t10 + w * t13)));
				p3 = t2 + w * (t5 + w * (t8 + w * (t11 + w * t14)));
				p = z * p1 - (tt - w * (p2 + y * p3));
				r += (tf + p);
				break;
			case 2:
				p1 = y * (u0 + y * (u1 + y * (u2 + y * (u3 + y * (u4 + y * u5)))));
				p2 = one + y * (v1 + y * (v2 + y * (v3 + y * (v4 + y * v5))));
				r += (-0.5 * y + p1 / p2);
			}
		} else if (ix < 0x40200000) { /* x < 8.0 */
			i = (int) x;
			t = zero;
			y = x - (double) i;
			p = y * (s0 + y * (s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y * s6))))));
			q = one + y * (r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * r6)))));
			r = half * y + p / q;
			z = one; /* lgamma(1+s) = log(s) + lgamma(s) */
			switch (i) {
			case 7:
				z *= (y + 6.0); /* FALLTHRU */
			case 6:
				z *= (y + 5.0); /* FALLTHRU */
			case 5:
				z *= (y + 4.0); /* FALLTHRU */
			case 4:
				z *= (y + 3.0); /* FALLTHRU */
			case 3:
				z *= (y + 2.0); /* FALLTHRU */
				r += Math.log(z);
				break;
			}
			/* 8.0 <= x < 2**58 */
		} else if (ix < 0x43900000) {
			t = Math.log(x);
			z = one / x;
			y = z * z;
			w = w0 + z * (w1 + y * (w2 + y * (w3 + y * (w4 + y * (w5 + y * w6)))));
			r = (x - half) * (t - one) + w;
		} else
			/* 2**58 <= x <= inf */
			r = x * (Math.log(x) - one);
		return r;
	}

	/**
	 * Calculates the log10 of the gamma function for x using the efficient
	 * FDLIBM implementation to avoid overflows and guarantees high accuracy
	 * even for large numbers.
	 *
	 * @param x
	 *            the x parameter
	 * @return the log10 of the gamma function at x.
	 */
	public static double log10Gamma(double x) {
		return lnToLog10(lnGamma(x));
	}

	public static double log10Factorial(int x) {
		if (x >= log10FactorialCache.length || x < 0)
			return log10Gamma(x + 1);
		else
			return log10FactorialCache[x];
	}

	/**
	 * Calculates the log10 of the binomial coefficient. Designed to prevent
	 * overflows even with very large numbers.
	 *
	 * @param n
	 *            total number of trials
	 * @param k
	 *            number of successes
	 * @return the log10 of the binomial coefficient
	 */
	public static double log10BinomialCoefficient(int n, int k) {
		return log10Factorial(n) - log10Factorial(k) - log10Factorial(n - k);
	}

	public static double log10BinomialProbability(int n, int k, double log10p) {
		double log10OneMinusP = Math.log10(1 - Math.pow(10, log10p));
		return log10BinomialCoefficient(n, k) + log10p * k + log10OneMinusP * (n - k);
	}

	/**
	 * Computes a binomial probability. This is computed using the formula
	 * <p/>
	 * B(k; n; p) = [ n! / ( k! (n - k)! ) ] (p^k)( (1-p)^k )
	 * <p/>
	 * where n is the number of trials, k is the number of successes, and p is
	 * the probability of success
	 *
	 * @param n
	 *            number of Bernoulli trials
	 * @param k
	 *            number of successes
	 * @param p
	 *            probability of success
	 * @return the binomial probability of the specified configuration. Computes
	 *         values down to about 1e-237.
	 */
	public static double binomialProbability(int n, int k, double p) {
		return Math.pow(10, log10BinomialProbability(n, k, Math.log10(p)));
	}

	public static int countOccurrences(char c, String s) {
		int count = 0;
		for (int i = 0; i < s.length(); i++) {
			count += s.charAt(i) == c ? 1 : 0;
		}
		return count;
	}

	public static int countOccurrences(final boolean element, final boolean[] array) {
		int count = 0;
		for (final boolean b : array) {
			if (element == b)
				count++;
		}

		return count;
	}

	/**
	 * calculate the Root Mean Square of an array of integers
	 *
	 * @param x
	 *            an int[] of numbers
	 * @return the RMS of the specified numbers.
	 */
	public static double rms(int[] x) {
		if (x.length == 0)
			return 0.0;

		double rms = 0.0;
		for (int i : x)
			rms += i * i;
		rms /= x.length;
		return Math.sqrt(rms);
	}

	/**
	 * A utility class that computes on the fly average and standard deviation
	 * for a stream of numbers. The number of observations does not have to be
	 * known in advance, and can be also very big (so that it could overflow any
	 * naive summation-based scheme or cause loss of precision). Instead, adding
	 * a new number <code>observed</code> to a sample with
	 * <code>add(observed)</code> immediately updates the instance of this
	 * object so that it contains correct mean and standard deviation for all
	 * the numbers seen so far. Source: Knuth, vol.2 (see also e.g.
	 * http://www.johndcook.com/standard_deviation.html for online reference).
	 */
	public static class RunningAverage {
		private double mean = 0.0;
		private double s = 0.0;
		private long obs_count = 0;

		public void add(double obs) {
			obs_count++;
			double oldMean = mean;
			mean += (obs - mean) / obs_count; // update mean
			s += (obs - oldMean) * (obs - mean);
		}

		public void addAll(Collection<Number> col) {
			for (Number o : col) {
				add(o.doubleValue());
			}
		}

		public double mean() {
			return mean;
		}

		public double stddev() {
			return Math.sqrt(s / (obs_count - 1));
		}

		public double var() {
			return s / (obs_count - 1);
		}

		public long observationCount() {
			return obs_count;
		}

		public RunningAverage clone() {
			RunningAverage ra = new RunningAverage();
			ra.mean = this.mean;
			ra.s = this.s;
			ra.obs_count = this.obs_count;
			return ra;
		}

		public void merge(RunningAverage other) {
			if (this.obs_count > 0 || other.obs_count > 0) { // if we have any
																// observations
																// at all
				this.mean = (this.mean * this.obs_count + other.mean * other.obs_count)
						/ (this.obs_count + other.obs_count);
				this.s += other.s;
			}
			this.obs_count += other.obs_count;
		}
	}

	public static double normalDistribution(final double mean, final double sd, final double x) {
		Utils.validateArg(sd >= 0, "sd: Standard deviation of normal must be >= 0");
		Utils.validateArg(wellFormedDouble(mean) && wellFormedDouble(sd) && wellFormedDouble(x),
				"mean, sd, or, x : Normal parameters must be well formatted (non-INF, non-NAN)");

		return (INV_SQRT_2_PI / sd) * Math.exp(-(x - mean) * (x - mean) / (2.0 * sd * sd));
	}

	public static double[] applyToArray(final double[] array, final DoubleUnaryOperator func) {
		Utils.nonNull(func, "function may not be null");
		Utils.nonNull(array, "array may not be null");
		final double[] result = new double[array.length];
		for (int m = 0; m < result.length; m++) {
			result[m] = func.applyAsDouble(array[m]);
		}
		return result;
	}

	public static double[] normalizeFromRealSpace(final double[] array) {
		if (array.length == 0)
			return array;

		final double sum = sum(array);
		Utils.validateArg(sum >= 0.0, () -> "Values in probability array sum to a negative number " + sum);
		return applyToArray(array, x -> x / sum);
	}

	public static void checkFinite(final double x) throws NotFiniteNumberException {
		if (Double.isInfinite(x) || Double.isNaN(x)) {
			throw new NotFiniteNumberException(x);
		}
	}
	
	public static double log10OneMinusX(final double x) {
	    if ( x == 1.0 )
	        return Double.NEGATIVE_INFINITY;
	    else if ( x == 0.0 )
	        return 0.0;
	    else {
	        final double d = Math.log10(1 / x - 1) + Math.log10(x);
	        return Double.isInfinite(d) || d > 0.0 ? 0.0 : d;
	    }
	}
	
	public static double log10BinomialProbability(final int n, final int k) {
	    return log10BinomialCoefficient(n, k) + (n * FAIR_BINOMIAL_PROB_LOG10_0_5);
	}

	public static int median(final int[] values) {
		Utils.nonNull(values);
		return (int) FastMath.round(new Median().evaluate(Arrays.stream(values).mapToDouble(n -> n).toArray()));
	}

	/**
	 * Compute the median of a list of numbers
	 *
	 * If values.length is even, this will be the middle value when the elements are sorted
	 * If values.length is odd then it will be the mean of the two values closest to the middle.
	 *
	 * @param values a list of numbers
	 * @return the median element of values
	 */
	public static <T extends Number & Comparable<T>> double median(final Collection<T> values) {
		Utils.nonEmpty(values, "cannot take the median of a collection with no values.");
		return new Median().evaluate(values.stream().mapToDouble(Number::doubleValue).toArray());
	}
}
