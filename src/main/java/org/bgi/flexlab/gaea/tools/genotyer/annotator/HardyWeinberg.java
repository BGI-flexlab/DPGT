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
package org.bgi.flexlab.gaea.tools.genotyer.annotator;

import htsjdk.tribble.util.popgen.HardyWeinbergCalculation;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.bgi.flexlab.gaea.data.structure.pileup.Mpileup;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.vcf.VariantDataTracker;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.WorkInProgressAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PerReadAlleleLikelihoodMap;
import org.bgi.flexlab.gaea.util.QualityUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Phred-scaled P value of genotype-based (using GT field) test for Hardy-Weinberg test for disequilibrium
 */
public class HardyWeinberg extends InfoFieldAnnotation implements WorkInProgressAnnotation {

    private static final int MIN_SAMPLES = 10;
    private static final int MIN_GENOTYPE_QUALITY = 10;
    private static final int MIN_LOG10_PERROR = MIN_GENOTYPE_QUALITY / 10;

    public Map<String, Object> annotate(final VariantDataTracker tracker,
                                        final ChromosomeInformationShare ref,
                                        final Mpileup mpileup,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {

        final GenotypesContext genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() < MIN_SAMPLES )
            return null;

        int refCount = 0;
        int hetCount = 0;
        int homCount = 0;
        for ( final Genotype g : genotypes ) {
            if ( g.isNoCall() )
                continue;

            // TODO - fix me:
            // Right now we just ignore genotypes that are not confident, but this throws off
            //  our HW ratios.  More analysis is needed to determine the right thing to do when
            //  the genotyper cannot decide whether a given sample is het or hom var.
            if ( g.getLog10PError() > MIN_LOG10_PERROR )
                continue;

            if ( g.isHomRef() )
                refCount++;
            else if ( g.isHet() )
                hetCount++;
            else
                homCount++;
        }

        if ( refCount + hetCount + homCount == 0)
            return null;

        double pvalue = HardyWeinbergCalculation.hwCalculate(refCount, hetCount, homCount);
        //System.out.println(refCount + " " + hetCount + " " + homCount + " " + pvalue);
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.1f", QualityUtils.phredScaleErrorRate(pvalue)));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("HW"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("HW", 1, VCFHeaderLineType.Float, "Phred-scaled p-value for Hardy-Weinberg violation")); }
}