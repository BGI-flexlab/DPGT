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
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.header.VCFConstants;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.vcf.VariantDataTracker;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.GenotypeAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PerReadAlleleLikelihoodMap;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DepthPerAlleleBySample extends GenotypeAnnotation implements StandardAnnotation {

    public void annotate(final VariantDataTracker tracker,
                         final ChromosomeInformationShare ref,
                         final Pileup pileup,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final PerReadAlleleLikelihoodMap alleleLikelihoodMap) {
        if ( g == null || !g.isCalled() || ( pileup == null && alleleLikelihoodMap == null) )
            return;

        if (alleleLikelihoodMap != null && !alleleLikelihoodMap.isEmpty())
            annotateWithLikelihoods(alleleLikelihoodMap, vc, gb);
        else if ( pileup != null && (vc.isSNP()))
            annotateWithPileup(pileup, vc, gb);
    }

    private void annotateWithPileup(final Pileup pileup, final VariantContext vc, final GenotypeBuilder gb) {
        //System.out.println("annotated with pileup");
        HashMap<Byte, Integer> alleleCounts = new HashMap<Byte, Integer>();
        for ( Allele allele : vc.getAlleles() ) {
            alleleCounts.put(allele.getBases()[0], 0);
            //System.out.println("genotype allele:" + allele.getBases()[0]);
        }

        for ( PileupReadInfo p : pileup.getTotalPileup() ) {
            if ( alleleCounts.containsKey(p.getByteBase()) ) {
                //System.out.println("pileup allele:" + p.getBase());
                alleleCounts.put(p.getByteBase(), alleleCounts.get(p.getByteBase())+1);
            }
        }

        // we need to add counts in the correct order
        int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference().getBases()[0]);
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
            counts[i+1] = alleleCounts.get(vc.getAlternateAllele(i).getBases()[0]);

        gb.AD(counts);
    }

    private void annotateWithLikelihoods(final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap, final VariantContext vc, final GenotypeBuilder gb) {
        final HashMap<Allele, Integer> alleleCounts = new HashMap<Allele, Integer>();

       //System.out.println("vc allele:");
        for ( final Allele allele : vc.getAlleles() ) {
            alleleCounts.put(allele, 0);
            //System.out.println(allele);
        }
        Map<Allele, Allele> alleles = getExtendedAllele(perReadAlleleLikelihoodMap.getRefAllele(), vc);
       /* System.out.println("recal vc allele:");
        for ( final Allele allele : alleles.keySet() ) {
            System.out.println(allele);
        }*/
        
        for (Map.Entry<AlignmentsBasic,Map<Allele,Double>> el : perReadAlleleLikelihoodMap.getLikelihoodReadMap().entrySet()) {
            final AlignmentsBasic read = el.getKey();
            final Allele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
            //System.out.println(read.getReadName() + "\t" + a.getBaseString());
            if (a.isNoCall())
                continue; // read is non-informative
            if (!alleles.keySet().contains(a))
                continue; // sanity check - shouldn't be needed
            Allele ab = alleles.get(a);
            alleleCounts.put(ab, alleleCounts.get(ab) + 1);
        }
        final int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference());
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
            counts[i+1] = alleleCounts.get( vc.getAlternateAllele(i) );

        gb.AD(counts);
    }
    
    private Map<Allele, Allele> getExtendedAllele(Allele ref, VariantContext vc) {
    	Map<Allele, Allele> alleles = new HashMap<Allele, Allele>();
    	Allele vcRef = vc.getReference();
    	if(vcRef == ref)
    		for ( Allele a : vc.getAlleles() ) {
    			alleles.put(a, a);
    		}
    	else {
    		//extend
    		if(vcRef.length() < ref.length()) {
    			byte[] extraBases = Arrays.copyOfRange(ref.getBases(), vcRef.length(), ref.length());
    			for ( Allele a : vc.getAlleles() ) {
                    if ( a.isReference() ) {
                    	alleles.put(ref, a);
                    }
                    else {
                        Allele extended = Allele.extend(a, extraBases);
                        alleles.put(extended, a);
                    }
    			}
    		}
    	}
    	
    	return alleles;
    }

    public List<String> getKeyNames() { return Arrays.asList(VCFConstants.GENOTYPE_ALLELE_DEPTHS); }

    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getFormatLine(getKeyNames().get(0)));
    }
}