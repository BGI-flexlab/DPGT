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
package org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator;

import org.bgi.flexlab.gaea.util.BaseUtils;
import org.bgi.flexlab.gaea.util.DiploidGenotype;


/**
 * Created by zhangyong on 2016/12/27.
 */
public class GenotypeData {

    /**
     * Constant static data: genotype zeros
     */
    public final static double[] genotypeZeros = new double[DiploidGenotype.values().length];

    /**
     * Constant static data: base zeros
     */
    public final static double[] baseZeros = new double[BaseUtils.BASES.length];

    static {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            genotypeZeros[g.ordinal()] = 0.0;
        }
        for ( byte base : BaseUtils.BASES ) {
            baseZeros[BaseUtils.simpleBaseToBaseIndex(base)] = 0.0;
        }
    }

    /**
     * The fundamental data arrays associated with a Genotype Likelihoods object
     */
    protected double[] log10Likelihoods = null;

    public GenotypeData() {
        setToZero();
    }

    protected void setToZero() {
        log10Likelihoods = genotypeZeros.clone();                 // likelihoods are all zeros
    }

    public void add(GenotypeData genotypeData) {
        for(int i = 0; i < log10Likelihoods.length; i++) {
            log10Likelihoods[i] += genotypeData.log10Likelihoods[i];
        }
    }

    public double[] getLog10Likelihoods() {
        return log10Likelihoods;
    }
}
