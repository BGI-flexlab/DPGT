package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.DoublePredicate;
import java.util.function.DoubleUnaryOperator;
import java.util.function.IntPredicate;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.util.MathUtils;

public class GvcfMathUtils extends MathUtils {

	private static final double LN_10 = Math.log(10);

	private static final double LOG1MEXP_THRESHOLD = Math.log(0.5);

	private final static double RAW_MIN_PHRED_SCALED_QUAL = Math.log10(Double.MIN_VALUE);

	private static final long GAEA_RANDOM_SEED = 47382911L;
	private static Random randomGenerator = new Random(GAEA_RANDOM_SEED);

	public static Random getRandomGenerator() {
		return randomGenerator;
	}

	public static void resetRandomGenerator() {
		randomGenerator.setSeed(GAEA_RANDOM_SEED);
	}

	public static void resetRandomGenerator(long seed) {
		randomGenerator.setSeed(seed);
	}

	private static final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(
			new Well19937c(GAEA_RANDOM_SEED));

	public static RandomDataGenerator getRandomDataGenerator() {
		return randomDataGenerator;
	}

	public static final double LOG10_OF_E = Math.log10(Math.E);

	public static final double LOG10_ONE_THIRD = -Math.log10(3.0);

	/**
	 * Compute Z=X-Y for two numeric vectors X and Y
	 *
	 * @param x
	 *            First vector
	 * @param y
	 *            Second vector
	 * @return Vector of same length as x and y so that z[k] = x[k]-y[k]
	 */
	public static int[] vectorDiff(final int[] x, final int[] y) {
		if (x.length != y.length)
			throw new UserException("BUG: Lengths of x and y must be the same");

		int[] result = new int[x.length];
		for (int k = 0; k < x.length; k++)
			result[k] = x[k] - y[k];

		return result;
	}

	public static double log10SumLog10(final double a, final double b) {
		return a > b ? a + Math.log10(1 + Math.pow(10.0, b - a)) : b + Math.log10(1 + Math.pow(10.0, a - b));
	}

	public static double log10OneMinusPow10(final double a) {
		if (a > 0)
			return Double.NaN;
		if (a == 0)
			return Double.NEGATIVE_INFINITY;
		final double b = a * LN_10;
		return log1mexp(b) / LN_10;
	}

	public static double log1mexp(final double a) {
		if (a > 0)
			return Double.NaN;
		if (a == 0)
			return Double.NEGATIVE_INFINITY;

		return (a < LOG1MEXP_THRESHOLD) ? Math.log1p(-Math.exp(a)) : Math.log(-Math.expm1(a));
	}

	/**
	 * Calculates the log10 of the multinomial coefficient. Designed to prevent
	 * overflows even with very large numbers.
	 *
	 * @param n
	 *            total number of trials
	 * @param k
	 *            array of any size with the number of successes for each
	 *            grouping (k1, k2, k3, ..., km)
	 * @return {@link Double#NaN NaN} if {@code a > 0}, otherwise the
	 *         corresponding value.
	 */
	public static double log10MultinomialCoefficient(final int n, final int[] k) {
		if (n < 0)
			throw new IllegalArgumentException("n: Must have non-negative number of trials");
		double denominator = 0.0;
		int sum = 0;
		for (int x : k) {
			if (x < 0)
				throw new IllegalArgumentException("x element of k: Must have non-negative observations of group");
			if (x > n)
				throw new IllegalArgumentException("x element of k, n: Group observations must be bounded by k");
			denominator += log10Factorial(x);
			sum += x;
		}
		if (sum != n)
			throw new IllegalArgumentException(
					"k and n: Sum of observations in multinomial must sum to total number of trials");
		return log10Factorial(n) - denominator;
	}

	public static class Log10Cache {
		/**
		 * Get the value of log10(n), expanding the cache as necessary
		 * 
		 * @param n
		 *            operand
		 * @return log10(n)
		 */
		public static double get(final int n) {
			if (n < 0)
				throw new UserException(String.format("Can't take the log of a negative number: %d", n));
			if (n >= cache.length)
				ensureCacheContains(Math.max(n + 10, 2 * cache.length));
			/*
			 * Array lookups are not atomic. It's possible that the reference to
			 * cache could be changed between the time the reference is loaded
			 * and the data is fetched from the correct offset. However, the
			 * value retrieved can't change, and it's guaranteed to be present
			 * in the old reference by the conditional above.
			 */
			return cache[n];
		}

		/**
		 * Ensures that the cache contains a value for n. After completion of
		 * ensureCacheContains(n), #get(n) is guaranteed to return without
		 * causing a cache expansion
		 * 
		 * @param n
		 *            desired value to be precomputed
		 */
		public static synchronized void ensureCacheContains(final int n) {
			if (n < cache.length)
				return;
			final double[] newCache = new double[n + 1];
			System.arraycopy(cache, 0, newCache, 0, cache.length);
			for (int i = cache.length; i < newCache.length; i++)
				newCache[i] = Math.log10(i);
			cache = newCache;
		}

		// initialize with the special case: log10(0) = NEGATIVE_INFINITY
		private static double[] cache = new double[] { Double.NEGATIVE_INFINITY };
	}

	private static class JacobianLogTable {

		public static final double MAX_TOLERANCE = 8.0;

		public static double get(final double difference) {
			if (cache == null)
				initialize();
			final int index = fastRound(difference * INV_STEP);
			return cache[index];
		}

		private static synchronized void initialize() {
			if (cache == null) {
				final int tableSize = (int) (MAX_TOLERANCE / TABLE_STEP) + 1;
				cache = new double[tableSize];
				for (int k = 0; k < cache.length; k++)
					cache[k] = Math.log10(1.0 + Math.pow(10.0, -((double) k) * TABLE_STEP));
			}
		}

		private static final double TABLE_STEP = 0.0001;
		private static final double INV_STEP = 1.0 / TABLE_STEP;
		private static double[] cache = null;
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
		if (diff >= MAX_JACOBIAN_TOLERANCE)
			return big;

		// OK, so |y-x| < tol: we use the following identity then:
		// we need to compute log10(10^x + 10^y)
		// By Jacobian logarithm identity, this is equal to
		// max(x,y) + log10(1+10^-abs(x-y))
		// we compute the second term as a table lookup with integer
		// quantization
		// we have pre-stored correction for 0,0.1,0.2,... 10.0
		final int ind = fastRound(diff * JACOBIAN_LOG_TABLE_INV_STEP); // hard
																		// rounding
		return big + jacobianLogTable[ind];
	}

	public static double approximateLog10SumLog10(final double[] vals) {
		return approximateLog10SumLog10(vals, vals.length);
	}

	public static double approximateLog10SumLog10(final double[] vals, final int endIndex) {

		final int maxElementIndex = MathUtils.maxElementIndex(vals, endIndex);
		double approxSum = vals[maxElementIndex];

		for (int i = 0; i < endIndex; i++) {
			if (i == maxElementIndex || vals[i] == Double.NEGATIVE_INFINITY)
				continue;

			final double diff = approxSum - vals[i];
			if (diff < JacobianLogTable.MAX_TOLERANCE) {
				// See notes from the 2-inout implementation below
				approxSum += JacobianLogTable.get(diff);
			}
		}

		return approxSum;
	}

	public static double phredScaleErrorRate(final double errorRate) {
		return phredScaleLog10ErrorRate(Math.log10(errorRate));
	}

	public static double phredScaleLog10ErrorRate(final double errorRateLog10) {
		if (!goodLog10Probability(errorRateLog10))
			throw new IllegalArgumentException("errorRateLog10 must be good probability but got " + errorRateLog10);
		// abs is necessary for edge base with errorRateLog10 = 0 producing -0.0
		// doubles
		return Math.abs(-10.0 * Math.max(errorRateLog10, RAW_MIN_PHRED_SCALED_QUAL));
	}

	public static double binomialCoefficient(final int n, final int k) {
		return Math.pow(10, log10BinomialCoefficient(n, k));
	}

	public static double[] applyToArrayInPlace(final double[] array, final DoubleUnaryOperator func) {
		JointCallingUtils.nonNull(array, "array may not be null");
		JointCallingUtils.nonNull(func, "function may not be null");
		for (int m = 0; m < array.length; m++) {
			array[m] = func.applyAsDouble(array[m]);
		}
		return array;
	}

	/**
	 * Test whether all elements of a double[] array satisfy a double -> boolean
	 * predicate
	 */
	public static boolean allMatch(final double[] array, final DoublePredicate pred) {
		JointCallingUtils.nonNull(array, "array may not be null");
		JointCallingUtils.nonNull(pred, "predicate may not be null");
		for (final double x : array) {
			if (!pred.test(x)) {
				return false;
			}
		}
		return true;
	}

	public static boolean allMatch(final int[] array, final IntPredicate pred) {
		JointCallingUtils.nonNull(array, "array may not be null");
		JointCallingUtils.nonNull(pred, "predicate may not be null");
		for (final int x : array) {
			if (!pred.test(x)) {
				return false;
			}
		}
		return true;
	}

	public static int[] sampleIndicesWithoutReplacement(final int n, final int k) {
		// No error checking : RandomDataGenetator.nextPermutation does it
		return GvcfMathUtils.getRandomDataGenerator().nextPermutation(n, k);
	}

	public static double log10(int i) {
		return Log10Cache.get(i);
	}

	public static double[] scaleLogSpaceArrayForNumericalStability(final double[] array) {
		final double maxValue = arrayMax(array);
		return applyToArrayInPlace(array, x -> x - maxValue);
	}

	public static double log10SumLog10(final double[] log10Values) {
		return log10SumLog10(log10Values, 0);
	}

	public static double log10SumLog10(final double[] log10Values, final int start) {
		return log10SumLog10(log10Values, start, log10Values.length);
	}

	public static double log10SumLog10(final double[] log10Values, final int start, final int finish) {
		if (start >= finish) {
			return Double.NEGATIVE_INFINITY;
		}
		final int maxElementIndex = maxElementIndex(log10Values, start, finish);
		final double maxValue = log10Values[maxElementIndex];
		if (maxValue == Double.NEGATIVE_INFINITY) {
			return maxValue;
		}
		double sum = 1.0;
		for (int i = start; i < finish; i++) {
			final double curVal = log10Values[i];
			if (i == maxElementIndex || curVal == Double.NEGATIVE_INFINITY) {
				continue;
			} else {
				final double scaled_val = curVal - maxValue;
				sum += Math.pow(10.0, scaled_val);
			}
		}
		if (Double.isNaN(sum) || sum == Double.POSITIVE_INFINITY) {
			throw new IllegalArgumentException("log10 p: Values must be non-infinite and non-NAN");
		}
		return maxValue + (sum != 1.0 ? Math.log10(sum) : 0.0);
	}

	public static double logToLog10(final double ln) {
		return ln * LOG10_OF_E;
	}

	public static double mean(final double[] in, final int start, final int stop) {
		return stop <= start ? Double.NaN : Arrays.stream(in, start, stop).average().getAsDouble();
	}

	public static double[] normalizeFromLog10ToLinearSpace(final double[] array) {
		return normalizeLog10(array, false, true);
	}

	public static double[] normalizeLog10(final double[] array, final boolean takeLog10OfOutput,
			final boolean inPlace) {
		final double log10Sum = log10SumLog10(array);
		final double[] result = inPlace ? applyToArrayInPlace(array, x -> x - log10Sum)
				: applyToArray(array, x -> x - log10Sum);
		return takeLog10OfOutput ? result : applyToArrayInPlace(result, x -> Math.pow(10.0, x));
	}
	public static <T extends Comparable<? super T>> T median(final List<T> array) {
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

	public static int arrayMin(final List<Integer> array) {
		if (array == null || array.isEmpty())
			throw new IllegalArgumentException("Array must be non-null and non-empty");
		int min = array.get(0);
		for (final int i : array)
			if (i < min)
				min = i;
		return min;
	}
}
