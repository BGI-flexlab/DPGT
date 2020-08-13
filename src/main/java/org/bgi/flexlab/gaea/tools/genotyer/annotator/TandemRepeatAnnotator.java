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


import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.bgi.flexlab.gaea.data.structure.pileup.Mpileup;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.vcf.VariantDataTracker;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.INDELGenotypeLikelihoodCalculator;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PerReadAlleleLikelihoodMap;
import org.bgi.flexlab.gaea.util.GaeaVariantContextUtils;
import org.bgi.flexlab.gaea.util.Pair;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class TandemRepeatAnnotator extends InfoFieldAnnotation implements StandardAnnotation {
    private static final String STR_PRESENT = "STR";
    private static final String REPEAT_UNIT_KEY = "RU";
    private static final String REPEATS_PER_ALLELE_KEY = "RPA";
    public Map<String, Object> annotate(final VariantDataTracker tracker,
                                        final ChromosomeInformationShare ref,
                                        final Mpileup mpileup,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( !vc.isIndel())
            return null;
        //System.err.println("variant:" + vc.getReference() + "," + vc.getAlternateAllele(0) +
        //        "\nforward seq:" + ref.getBaseSequence(mpileup.getPosition(), Math.min(ref.getLength(), mpileup.getPosition() + INDELGenotypeLikelihoodCalculator.REF_WIN_Extend)));
        Pair<List<Integer>,byte[]> result = GaeaVariantContextUtils.getNumTandemRepeatUnits(vc,
                ref.getGA4GHBaseBytes(mpileup.getPosition(), Math.min(ref.getLength() - 1, mpileup.getPosition() + INDELGenotypeLikelihoodCalculator.REF_WIN_Extend)));
        if (result == null)
            return null;

        byte[] repeatUnit = result.second;
        List<Integer> numUnits = result.first;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(STR_PRESENT,true);
        map.put(REPEAT_UNIT_KEY,new String(repeatUnit));
        map.put(REPEATS_PER_ALLELE_KEY, numUnits);

        return map;
    }

    protected static final String[] keyNames = {STR_PRESENT, REPEAT_UNIT_KEY,REPEATS_PER_ALLELE_KEY };
    protected static final VCFInfoHeaderLine[] descriptions = {
            new VCFInfoHeaderLine(STR_PRESENT, 0, VCFHeaderLineType.Flag, "Variant is a short tandem repeat"),
            new VCFInfoHeaderLine(REPEAT_UNIT_KEY, 1, VCFHeaderLineType.String, "Tandem repeat unit (bases)"),
            new VCFInfoHeaderLine(REPEATS_PER_ALLELE_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Number of times tandem repeat unit is repeated, for each allele (including reference)") };

    public List<String> getKeyNames() {
        return Arrays.asList(keyNames);
    }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(descriptions); }

}
