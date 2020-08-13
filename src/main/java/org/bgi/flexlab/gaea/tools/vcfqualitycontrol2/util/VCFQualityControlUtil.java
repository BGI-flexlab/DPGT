package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util;

import java.util.function.DoubleUnaryOperator;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.util.Utils;

public class VCFQualityControlUtil {
	private static final double SQUARE_ROOT_OF_TWO_TIMES_PI = Math.sqrt(2.0 * Math.PI);
	private static final double NATURAL_LOG_OF_TEN = Math.log(10.0);
	
	public static double[] normalizeLog10DeleteMePlease(final double[] array, final boolean takeLog10OfOutput) {
        final double maxValue = arrayMax(array);
        final double[] normalized = applyToArray(array, x -> Math.pow(10.0, x - maxValue));
        final double sum = sum(normalized);
        if (!takeLog10OfOutput) {
            return applyToArrayInPlace(normalized, x -> x/sum);
        } else {
            final double log10Sum = Math.log10(sum);
            return applyToArrayInPlace(array, x -> x - maxValue - log10Sum);
        }
    }
	
	public static double arrayMax(final double[] array) {
        return array[maxElementIndex(array)];
    }
	
	public static int maxElementIndex(final double[] array) {
        return maxElementIndex(array, array.length);
    }

    public static int maxElementIndex(final double[] array, final int start, final int endIndex) {
        Utils.nonNull(array, "array may not be null");
        if(array.length < 0)
        	throw new UserException("array may not be empty");
        if(start >= endIndex)
        	throw new UserException("Start cannot be after end.");

        int maxI = start;
        for (int i = (start+1); i < endIndex; i++) {
            if (array[i] > array[maxI])
                maxI = i;
        }
        return maxI;
    }

    public static int maxElementIndex(final double[] array, final int endIndex) {
        return maxElementIndex(array, 0, endIndex);
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
    
    public static double sum(final double[] values) {
        double s = 0.0;
        for (double v : values)
            s += v;
        return s;
    }
    
    public static double[] applyToArrayInPlace(final double[] array, final DoubleUnaryOperator func) {
        Utils.nonNull(array, "array may not be null");
        Utils.nonNull(func, "function may not be null");
        for (int m = 0; m < array.length; m++) {
            array[m] = func.applyAsDouble(array[m]);
        }
        return array;
    }
    
    public static double normalDistributionLog10(final double mean, final double sd, final double x) {
    	if(sd < 0)
    		throw new UserException("sd: Standard deviation of normal must be > 0");
        if ( ! wellFormedDouble(mean) || ! wellFormedDouble(sd) || ! wellFormedDouble(x) )
            throw new IllegalArgumentException("mean, sd, or, x : Normal parameters must be well formatted (non-INF, non-NAN)");
        final double a = -1.0 * Math.log10(sd * SQUARE_ROOT_OF_TWO_TIMES_PI);
        final double b = -1.0 * (square(x - mean) / (2.0 * square(sd))) / NATURAL_LOG_OF_TEN;
        return a + b;
    }
    
    public static boolean wellFormedDouble(final double val) {
        return !Double.isInfinite(val) && !Double.isNaN(val);
    }
    
    public static double square(final double x) {
        return x * x;
    }
}
