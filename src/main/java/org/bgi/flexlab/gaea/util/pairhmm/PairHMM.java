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
package org.bgi.flexlab.gaea.util.pairhmm;

//import com.google.java.contract.Ensures;
//import com.google.java.contract.Requires;


/**
 * from GATK
 */
public abstract class PairHMM {
    protected static final Byte MAX_CACHED_QUAL = Byte.MAX_VALUE;
    protected static final byte DEFAULT_GOP = (byte) 45;
    protected static final byte DEFAULT_GCP = (byte) 10;

    public enum HMM_IMPLEMENTATION {
        /* Very slow implementation which uses very accurate log10 sum functions. Only meant to be used as a reference test implementation */
        EXACT,
        /* PairHMM as implemented for the UnifiedGenotyper. Uses log10 sum functions accurate to only 1E-4 */
        ORIGINAL,
        /* Optimized version of the PairHMM which caches per-read computations */
        CACHING,
        /* Optimized version of the PairHMM which caches per-read computations and operations in real space to avoid costly sums of log10'ed likelihoods */
        LOGLESS_CACHING
    }

    protected double[][] matchMetricArray = null;
    protected double[][] XMetricArray = null;
    protected double[][] YMetricArray = null;

    public abstract void initialize( final int READ_MAX_LENGTH, final int HAPLOTYPE_MAX_LENGTH );

 //   @Requires({"readBases.length == readQuals.length", "readBases.length == insertionGOP.length", "readBases.length == deletionGOP.length",
 //              "readBases.length == overallGCP.length", "matchMetricArray!=null", "XMetricArray!=null", "YMetricArray!=null"})
 //   @Ensures({"!Double.isInfinite(result)", "!Double.isNaN(result)"}) // Result should be a proper log10 likelihood
    public abstract double computeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                                     final byte[] readBases,
                                                                     final byte[] readQuals,
                                                                     final byte[] insertionGOP,
                                                                     final byte[] deletionGOP,
                                                                     final byte[] overallGCP,
                                                                     final int hapStartIndex,
                                                                     final boolean recacheReadValues );
}
