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
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.vcf.VariantDataTracker;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.ExperimentalAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.GenotypeAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PerReadAlleleLikelihoodMap;
import org.bgi.flexlab.gaea.util.MathUtils;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * The allele balance (fraction of ref bases over ref + alt bases) separately for each bialleleic het-called sample
 */
public class AlleleBalanceBySample extends GenotypeAnnotation implements ExperimentalAnnotation {

    public void annotate(final VariantDataTracker tracker,
                         final ChromosomeInformationShare ref,
                         final Pileup pileup,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final PerReadAlleleLikelihoodMap alleleLikelihoodMap){
        if ( pileup == null )
            return;

        Double ratio = annotateSNP(pileup, vc, g);
        if (ratio == null)
            return;

        gb.attribute(getKeyNames().get(0), Double.valueOf(String.format("%.2f", ratio.doubleValue())));
    }

    private Double annotateSNP(Pileup pileup, VariantContext vc, Genotype g) {
        double ratio = -1;

        if ( !vc.isSNP() )
            return null;

        if ( !vc.isBiallelic() )
            return null;

        if ( g == null || !g.isCalled() )
            return null;

        if (!g.isHet())
            return null;

        Collection<Allele> altAlleles = vc.getAlternateAlleles();
        if ( altAlleles.size() == 0 )
            return null;

        final String bases = new String(pileup.getBases());
        if ( bases.length() == 0 )
            return null;
        char refChr = vc.getReference().toString().charAt(0);
        char altChr = vc.getAlternateAllele(0).toString().charAt(0);

        int refCount = MathUtils.countOccurrences(refChr, bases);
        int altCount = MathUtils.countOccurrences(altChr, bases);

        // sanity check
        if ( refCount + altCount == 0 )
            return null;

        ratio = ((double)refCount / (double)(refCount + altCount));
        return ratio;
    }

    public List<String> getKeyNames() { return Arrays.asList("AB"); }

    public List<VCFFormatHeaderLine> getDescriptions() { return Arrays.asList(new VCFFormatHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Float, "Allele balance for each het genotype")); }
}