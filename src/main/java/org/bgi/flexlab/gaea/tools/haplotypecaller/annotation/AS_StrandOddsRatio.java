package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.AlleleSpecificAnnotationData;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public final class AS_StrandOddsRatio extends AS_StrandBiasTest implements AS_StandardAnnotation {

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GaeaVCFConstants.AS_STRAND_ODDS_RATIO_KEY);
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromLikelihoods(final ReadLikelihoods<Allele> likelihoods,
                                                                     final VariantContext vc){
        // either SNP with no alignment context, or indels: per-read likelihood map needed
        final int[][] table = getContingencyTable(likelihoods, vc, MIN_COUNT);
        final double ratio = StrandOddsRatio.calculateSOR(table);
        return Collections.singletonMap(getKeyNames().get(0), StrandOddsRatio.formattedValue(ratio));
    }

    @Override
    protected Map<Allele,Double> calculateReducedData(AlleleSpecificAnnotationData<List<Integer>> combinedData) {
        final Map<Allele,Double> annotationMap = new HashMap<>();
        final Map<Allele, List<Integer>> perAlleleData = combinedData.getAttributeMap();
        final List<Integer> refStrandCounts = perAlleleData.get(combinedData.getRefAllele());
        for (final Allele a : perAlleleData.keySet()) {
            List<Integer> altStrandCounts = perAlleleData.get(a);
            int[][] refAltTable = new int[][] {new int[]{refStrandCounts.get(FORWARD),refStrandCounts.get(REVERSE)},
                    new int[]{altStrandCounts.get(FORWARD),altStrandCounts.get(REVERSE)}};
            annotationMap.put(a,StrandOddsRatio.calculateSOR(refAltTable));
        }
        return annotationMap;
    }

}

