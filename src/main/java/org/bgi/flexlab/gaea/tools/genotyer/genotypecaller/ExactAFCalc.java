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

package org.bgi.flexlab.gaea.tools.genotyer.genotypecaller;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.bgi.flexlab.gaea.util.GaeaVariantContextUtils;
import org.bgi.flexlab.gaea.util.MathUtils;

import java.util.ArrayList;

abstract class ExactAFCalc extends AlleleFrequencyCalculator {
    protected static final int HOM_REF_INDEX = 0;  // AA likelihoods are always first

    protected ExactAFCalc(final int nSamples, int maxAltAlleles, final int ploidy) {
        super(nSamples, maxAltAlleles, ploidy);
    }

    /**
     * Wrapper class that compares two likelihoods associated with two alleles
     */
    protected static final class LikelihoodSum implements Comparable<LikelihoodSum> {
        public double sum = 0.0;
        public Allele allele;

        public LikelihoodSum(Allele allele) { this.allele = allele; }

        public int compareTo(LikelihoodSum other) {
            final double diff = sum - other.sum;
            return ( diff < 0.0 ) ? 1 : (diff > 0.0 ) ? -1 : 0;
        }
    }

    /**
     * Unpack GenotypesContext into arraylist of doubel values
     * @param GLs            Input genotype context
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    protected static ArrayList<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy) {
        ArrayList<double[]> genotypeLikelihoods = new ArrayList<double[]>(GLs.size() + 1);
        if ( includeDummy ) genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ) {
        	if ( sample.hasLikelihoods() ) {
                double[] gls = sample.getLikelihoods().getAsVector();
                if ( MathUtils.sum(gls) < GaeaVariantContextUtils.SUM_GL_THRESH_NOCALL )
                    genotypeLikelihoods.add(gls);
            }
        }

        return genotypeLikelihoods;
    }
}