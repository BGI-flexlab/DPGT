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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PerReadAlleleLikelihoodMap;

import java.util.Arrays;
import java.util.List;
import java.util.Map;


/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 6/28/12
 */

/**
 * The u-based z-approximation from the Mann-Whitney Rank Sum Test for reads with clipped bases (reads with ref bases vs. those with the alternate allele)
 * Note that the clipping rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.
 */
public class ClippingRankSumTest extends RankSumTest {

    public List<String> getKeyNames() { return Arrays.asList("ClippingRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("ClippingRankSum", 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases")); }


    protected void fillQualsFromPileup(final List<Allele> allAlleles,
                                       final int refLoc,
                                       final Pileup pileup,
                                       final PerReadAlleleLikelihoodMap likelihoodMap, final List<Double> refQuals, final List<Double> altQuals) {
        // todo - only support non-pileup case for now, e.g. active-region based version
        if (pileup != null || likelihoodMap == null)
            return;

        for (Map.Entry<AlignmentsBasic, Map<Allele,Double>> el : likelihoodMap.getLikelihoodReadMap().entrySet()) {

            final Allele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
            if (a.isNoCall())
                continue; // read is non-informative
            if (a.isReference())
                refQuals.add((double)el.getKey().getNumHardClippedBases());
            else if (allAlleles.contains(a))
                altQuals.add((double)el.getKey().getNumHardClippedBases());

        }
    }

 }
