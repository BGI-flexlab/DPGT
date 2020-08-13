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
package org.bgi.flexlab.gaea.data.structure.variant;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Created by zhangyong on 2017/2/22.
 */
public class VariantCallContext extends VariantContext{
    // Was the site called confidently, either reference or variant?
    public boolean confidentlyCalled = false;

    // Should this site be emitted?
    public boolean shouldEmit = true;

    public VariantCallContext(VariantContext vc, boolean confidentlyCalledP) {
        super(vc);
        this.confidentlyCalled = confidentlyCalledP;
    }

    public VariantCallContext(VariantContext vc, boolean confidentlyCalledP, boolean shouldEmit) {
        super(vc);
        this.confidentlyCalled = confidentlyCalledP;
        this.shouldEmit = shouldEmit;
    }

    /* these methods are only implemented for GENOTYPE_GIVEN_ALLELES MODE */
    //todo -- expand these methods to all modes

    /**
     *
     * @param callConfidenceThreshold the Unified Argument Collection STANDARD_CONFIDENCE_FOR_CALLING
     * @return true if call was confidently ref
     */
    public boolean isCalledRef(double callConfidenceThreshold) {
        return (confidentlyCalled && (getPhredScaledQual() < callConfidenceThreshold));
    }

    /**
     *
     * @param callConfidenceThreshold the Unified Argument Collection STANDARD_CONFIDENCE_FOR_CALLING
     * @return true if call was confidently alt
     */
    public boolean isCalledAlt(double callConfidenceThreshold) {
        return (confidentlyCalled && (getPhredScaledQual() > callConfidenceThreshold));
    }
}
