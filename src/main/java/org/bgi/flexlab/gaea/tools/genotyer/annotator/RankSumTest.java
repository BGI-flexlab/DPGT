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
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.bgi.flexlab.gaea.data.structure.pileup.Mpileup;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.vcf.VariantDataTracker;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PerReadAlleleLikelihoodMap;
import org.bgi.flexlab.gaea.util.MannWhitneyU;
import org.bgi.flexlab.gaea.util.Pair;
import org.bgi.flexlab.gaea.util.QualityUtils;

import java.util.*;


/**
 * Abstract root for all RankSum based annotations
 */
public abstract class RankSumTest extends InfoFieldAnnotation implements ActiveRegionBasedAnnotation {
    //static final boolean DEBUG = false;
    private boolean useDithering = true;

    public Map<String, Object> annotate(final VariantDataTracker tracker,
                                        final ChromosomeInformationShare ref,
                                        final Mpileup mpileup,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        // either stratifiedContexts or  stratifiedPerReadAlleleLikelihoodMap has to be non-null
    	
        final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || genotypes.size() == 0)
            return null;

        final ArrayList<Double> refQuals = new ArrayList<Double>();
        final ArrayList<Double> altQuals = new ArrayList<Double>();

        for ( final Genotype genotype : genotypes.iterateInSampleNameOrder() ) {
            PerReadAlleleLikelihoodMap indelLikelihoodMap = null;

            if (stratifiedPerReadAlleleLikelihoodMap != null )
                indelLikelihoodMap = stratifiedPerReadAlleleLikelihoodMap.get(genotype.getSampleName());

            if (indelLikelihoodMap != null && indelLikelihoodMap.isEmpty())
                indelLikelihoodMap = null;
            // treat an empty likelihood map as a null reference - will simplify contract with fillQualsFromPileup
            if (indelLikelihoodMap == null && mpileup == null)
                continue;

            Pileup pileup = mpileup.getCurrentPosPileup().get(genotype.getSampleName());
            fillQualsFromPileup(vc.getAlleles(), vc.getStart(), pileup, indelLikelihoodMap, refQuals, altQuals );
        }
        
        if (refQuals.isEmpty() && altQuals.isEmpty())
            return null;

        final MannWhitneyU mannWhitneyU = new MannWhitneyU(useDithering);
        //System.err.println("useDithering:"+useDithering);
        for (final Double qual : altQuals) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET1);
            //System.err.println("add alt:"+qual);
        }
        for (final Double qual : refQuals) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET2);
            //System.err.println("add ref:"+qual);
        }

       /* System.out.format("%s, REF QUALS:", this.getClass().getName());
        for (final Double qual : refQuals)
            System.out.format("%4.1f ", qual);
        System.out.println();
        System.out.format("%s, ALT QUALS:", this.getClass().getName());
        for (final Double qual : altQuals)
            System.out.format("%4.1f ", qual);
        System.out.println();*/
        
        // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
        final Pair<Double, Double> testResults = mannWhitneyU.runOneSidedTest(MannWhitneyU.USet.SET1);
        //System.err.println(getKeyNames().get(0));
        //System.err.println("testResult:"+testResults.getFirst()+"\t"+testResults.getSecond());
        final Map<String, Object> map = new HashMap<String, Object>();
        if (!Double.isNaN(testResults.first))
            map.put(getKeyNames().get(0), String.format("%.3f", testResults.first));
        return map;
    }

     protected abstract void fillQualsFromPileup(final List<Allele> alleles,
                                                final int refLoc,
                                                final Pileup pileup,
                                                final PerReadAlleleLikelihoodMap alleleLikelihoodMap,
                                                final List<Double> refQuals,
                                                final List<Double> altQuals);

    /**
     * Can the base in this pileup element be used in comparative tests between ref / alt bases?
     *
     * Note that this function by default does not allow deletion pileup elements
     *
     * @param p the pileup element to consider
     * @return true if this base is part of a meaningful read for comparison, false otherwise
     */
    public static boolean isUsableBase(final PileupReadInfo p) {
        return isUsableBase(p, false);
    }

    /**
     * Can the base in this pileup element be used in comparative tests between ref / alt bases?
     *
     * @param p the pileup element to consider
     * @param allowDeletions if true, allow p to be a deletion base
     * @return true if this base is part of a meaningful read for comparison, false otherwise
     */
    public static boolean isUsableBase(final PileupReadInfo p, final boolean allowDeletions) {
        return !(p.isInsertionAtBeginningOfRead() ||
                 (! allowDeletions && p.isDeletionBase()) ||
                 p.getReadInfo().getMappingQual() == 0 ||
                p.getReadInfo().getMappingQual() == QualityUtils.MAPPING_QUALITY_UNAVAILABLE ||
                 ((int) p.getBaseQuality()) < QualityUtils.MINIMUM_USABLE_QUALITY_SCORE); // need the unBAQed quality score here
    }

    /**
     * Initialize the rank sum test annotation using walker and engine information. Right now this checks to see if
     * engine randomization is turned off, and if so does not dither.
     * @param headerLines
     */
    public void initialize ( Set<VCFHeaderLine> headerLines ) {
       // useDithering = ! toolkit.getArguments().disableRandomization;
    	useDithering=true;
    }
    
   
}