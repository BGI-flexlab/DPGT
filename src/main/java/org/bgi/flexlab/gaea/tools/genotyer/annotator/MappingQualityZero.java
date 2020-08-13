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
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.pileup.Mpileup;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.vcf.VariantDataTracker;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PerReadAlleleLikelihoodMap;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Total count across all samples of mapping quality zero reads
 */
public class MappingQualityZero extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

    public Map<String, Object> annotate(final VariantDataTracker tracker,
                                        final ChromosomeInformationShare ref,
                                        final Mpileup mpileup,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ((vc.isSNP() || !vc.isVariant()) && mpileup != null)
            return annotatePileup(ref, mpileup, vc);
        else if (stratifiedPerReadAlleleLikelihoodMap != null && vc.isVariant())
            return annotateWithLikelihoods(stratifiedPerReadAlleleLikelihoodMap, vc);
        else
            return null;
    }

    private Map<String, Object> annotatePileup(final ChromosomeInformationShare ref,
                                               final Mpileup mpileup,
                                               final VariantContext vc) {
        if ( mpileup.getSize() == 0 )
            return null;

        int mq0 = 0;
        for ( String sample : mpileup.getCurrentPosPileup().keySet() ) {
            Pileup pileup = mpileup.getCurrentPosPileup().get(sample);

            for (PileupReadInfo p : pileup.getTotalPileup() ) {
                if ( p.getMappingQuality() == 0 )
                    mq0++;
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", mq0));
        return map;
    }

    private Map<String, Object> annotateWithLikelihoods(final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                        final VariantContext vc) {
        if (stratifiedPerReadAlleleLikelihoodMap == null)
            return null;

        int mq0 = 0;
        for ( PerReadAlleleLikelihoodMap likelihoodMap : stratifiedPerReadAlleleLikelihoodMap.values() ) {
            for (AlignmentsBasic read : likelihoodMap.getLikelihoodReadMap().keySet()) {

                 if (read.getMappingQual() == 0 )
                    mq0++;
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", mq0));
        return map;
    }


    public List<String> getKeyNames() { return Arrays.asList(VCFConstants.MAPPING_QUALITY_ZERO_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}