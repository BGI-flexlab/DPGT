package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.AlleleSpecificAnnotationData;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.QualityUtils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class AS_FisherStrand extends AS_StrandBiasTest implements AS_StandardAnnotation {

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GaeaVCFConstants.AS_FISHER_STRAND_KEY);
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromLikelihoods(final ReadLikelihoods<Allele> likelihoods,
                                                                     final VariantContext vc) {
        // either SNP with no alignment context, or indels: per-read likelihood map needed
        final int[][] table = StrandBiasTest.getContingencyTable(likelihoods, vc, MIN_COUNT);
        return table == null ? null : annotationForOneTable(FisherStrand.pValueForContingencyTable(table));
    }

    /**
     * Returns an annotation result given a pValue
     *
     * @return a hash map from FS -> phred-scaled pValue
     */
    private Map<String, Object> annotationForOneTable(final double pValue) {
        return Collections.singletonMap(getKeyNames().get(0), FisherStrand.makeValueObjectForAnnotation(pValue));
    }

    @Override
    protected Map<Allele,Double> calculateReducedData(AlleleSpecificAnnotationData<List<Integer>> combinedData) {
        final Map<Allele,Double> annotationMap = new HashMap<>();
        final Map<Allele,List<Integer>> perAlleleData = combinedData.getAttributeMap();
        final List<Integer> refStrandCounts = perAlleleData.get(combinedData.getRefAllele());
        for (final Allele a : perAlleleData.keySet()) {
            if(!a.equals(combinedData.getRefAllele(),true)) {
                final List<Integer> altStrandCounts = combinedData.getAttribute(a);
                final int[][] refAltTable = new int[][]{new int[]{refStrandCounts.get(0), refStrandCounts.get(1)}, new int[]{altStrandCounts.get(0), altStrandCounts.get(1)}};
                annotationMap.put(a, QualityUtils.phredScaleErrorRate(Math.max(FisherStrand.pValueForContingencyTable(refAltTable), MIN_PVALUE)));
            }
        }
        return annotationMap;
    }

}

