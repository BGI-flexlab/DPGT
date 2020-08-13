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

public class QualityUtils {
	public final static byte MAXIMUM_QUALITY_SCORE = 93;
	public final static byte MAXIMUM_USABLE_QUALITY_SCORE = 93;
	public final static byte MINIMUM_USABLE_QUALITY_SCORE = 6;
	public final static int MAPPING_QUALITY_UNAVAILABLE = 255;

	private static double qualityToErrorProbabilityCache[] = new double[256];
	public static double[] QUALITY_PROB = new double[MAXIMUM_USABLE_QUALITY_SCORE + 1];
	public static double[] QUALITY_PROB_LOG10 = new double[MAXIMUM_USABLE_QUALITY_SCORE + 1];
	public static double[] MINUS_QUALITY_PROB_LOG10 = new double[MAXIMUM_USABLE_QUALITY_SCORE + 1];
	
	private static final double[] qualToErrorProbCache = new double[MAXIMUM_QUALITY_SCORE + 1];
    private static final double[] qualToProbLog10Cache = new double[MAXIMUM_QUALITY_SCORE + 1];


    static {
        for (int i = 0; i <= MAXIMUM_QUALITY_SCORE; i++) {
            qualToErrorProbCache[i] = qualToErrorProb((double) i);
            qualToProbLog10Cache[i] = Math.log10(1.0 - qualToErrorProbCache[i]);
        }
    }
    
    public static double qualToErrorProb(final double qual) {
        Utils.validateArg( qual >= 0.0, () -> "qual must be >= 0.0 but got " + qual);
        return Math.pow(10.0, qual / -10.0);
    }
	
	private static final double MIN_LOG10_SCALED_QUAL = Math.log10(Double.MIN_VALUE);

	static {
		for(byte quality = 0; quality < MAXIMUM_USABLE_QUALITY_SCORE + 1; quality++) {
			double errorProbability = qualityToErrorProbability(quality);
			QUALITY_PROB[quality] = errorProbability;
			QUALITY_PROB_LOG10[quality] = quality / -10.0;
			MINUS_QUALITY_PROB_LOG10[quality] = Math.log10(1.0 - errorProbability);
		}
	}

	static {
		for (int i = 0; i < 256; i++)
			qualityToErrorProbabilityCache[i] = qualityToErrorProbability(i);
	}

	static public double qualToProbLog10(byte qual) {
		return MINUS_QUALITY_PROB_LOG10[(int)qual & 0xff];
	}

	static public double qualToErrorProbLog10(byte qual) {
		return QUALITY_PROB_LOG10[(int)qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
	}
	
	public static double qualToErrorProbLog10(final double qual) {
	    Utils.validateArg( qual >= 0.0, () -> "qual must be >= 0.0 but got " + qual);
	    return qual / -10.0;
	}

	public static double qualityToErrorProbability(final double qual) {
		return Math.pow(10.0, qual / -10.0);
	}
	
	static public double qualityToProbability(double qual) {
        return 1.0 - qualityToErrorProbability(qual);
    }

	public static byte probabilityToQuality(double prob, double eps) {
		double lp = Math.round(-10.0 * Math.log10(1.0 - prob + eps));
		return boundQuality((int) lp);
	}

	public static double phredScaleErrorRate(double errorRate) {
		return Math.abs(-10.0 * Math.log10(errorRate));
	}

	public static byte boundQuality(int qual) {
		return boundQuality(qual, MAXIMUM_USABLE_QUALITY_SCORE);
	}

	public static byte boundQuality(int qual, byte maxQual) {
		return (byte) Math.max(Math.min(qual, maxQual), 1);
	}

    /**
     * Convert a log10 probability of being right to a phred-scaled quality score of being wrong as a double
     *
     * This is a very generic method, that simply computes a phred-scaled double quality
     * score given an error rate.  It has the same precision as a normal double operation
     */
    public static double phredScaleLog10CorrectRate(final double trueRateLog10) {
        return phredScaleCorrectRate(Math.pow(10.0, trueRateLog10));
    }
    
    public static double phredScaleCorrectRate(final double trueRate) {
        return phredScaleLog10ErrorRate(MathUtils.log10OneMinusX(trueRate));
    }
    
    public static double phredScaleLog10ErrorRate(final double errorRateLog10) {
        if ( ! MathUtils.goodLog10Probability(errorRateLog10) ) throw new
                IllegalArgumentException("errorRateLog10 must be good probability but got " + errorRateLog10);
        // abs is necessary for edge base with errorRateLog10 = 0 producing -0.0 doubles
        return Math.abs(-10.0 * Math.max(errorRateLog10, MIN_LOG10_SCALED_QUAL));
    }
    
    public static double qualToErrorProb(final byte qual) {
        return qualToErrorProbCache[(int)qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }
    
    public static double qualToProb(final byte qual) {
        return 1.0 - qualToErrorProb(qual);
    }
    
    public static double qualToProb(final double qual) {
        Utils.validateArg( qual >= 0.0, () -> "qual must be >= 0.0 but got " + qual);
        return 1.0 - qualToErrorProb(qual);
    }
}
