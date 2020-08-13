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
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PerReadAlleleLikelihoodMap;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * The u-based z-approximation from the Mann-Whitney Rank Sum Test for base qualities (ref bases vs. bases of the alternate allele).
 * Note that the base quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.
 */
public class BaseQualityRankSumTest extends RankSumTest implements StandardAnnotation {
    public List<String> getKeyNames() { return Arrays.asList("BaseQRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("BaseQRankSum", 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities")); }

    protected void fillQualsFromPileup(final List<Allele> allAlleles, final int refLoc,
                                       final Pileup pileup,
                                       final PerReadAlleleLikelihoodMap alleleLikelihoodMap,
                                       final List<Double> refQuals, final List<Double> altQuals){

        if (alleleLikelihoodMap == null) {
            // use fast SNP-based version if we don't have per-read allele likelihoods
            for ( final PileupReadInfo p : pileup.getTotalPileup() ) {
                if ( isUsableBase(p) ) {
                    if ( allAlleles.get(0).equals(Allele.create((byte)p.getBase(),true)) ) {
                        refQuals.add((double)p.getBaseQuality());
                    } else if ( allAlleles.contains(Allele.create((byte)p.getBase()))) {
                        altQuals.add((double)p.getBaseQuality());
                    }
                }
            }
            
            
            return;
        }

        for (Map<Allele,Double> el : alleleLikelihoodMap.getLikelihoodMapValues()) {
            final Allele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el);
            if (a.isNoCall())
                continue; // read is non-informative
            if (a.isReference())
                refQuals.add(-10.0*(double)el.get(a));
            else if (allAlleles.contains(a))
                altQuals.add(-10.0*(double)el.get(a));


        }
       
    }


}