package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

import static java.lang.Math.max;
import static java.lang.Math.min;

public final class StrandOddsRatio extends StrandBiasTest implements StandardAnnotation {

    private static final double PSEUDOCOUNT = 1.0;
    private static final int MIN_COUNT = 0;

    @Override
    protected Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes){
        final int[][] tableFromPerSampleAnnotations = getTableFromSamples(genotypes, MIN_COUNT);
        return tableFromPerSampleAnnotations != null ? annotationForOneTable(calculateSOR(tableFromPerSampleAnnotations)) : null;
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromLikelihoods(final ReadLikelihoods<Allele> likelihoods, final VariantContext vc){
        final int[][] table = getContingencyTable(likelihoods, vc, MIN_COUNT);
        return annotationForOneTable(calculateSOR(table));
    }

    /**
     * Computes the SOR value of a table after augmentation. Based on the symmetric odds ratio but modified to take on
     * low values when the reference +/- read count ratio is skewed but the alt count ratio is not.  Natural log is taken
     * to keep values within roughly the same range as other annotations.
     *
     * Adding pseudocounts avoids division by zero.
     *
     * @param table The table before adding pseudocounts
     * @return the SOR annotation value
     */
    public static double calculateSOR(final int[][] table) {
        final double t00 = table[0][0] + PSEUDOCOUNT;
        final double t01 = table[0][1] + PSEUDOCOUNT;
        final double t11 = table[1][1] + PSEUDOCOUNT;
        final double t10 = table[1][0] + PSEUDOCOUNT;

        final double ratio = (t00 / t01) * (t11 / t10) + (t01 / t00) * (t10 / t11);

        final double refRatio = min(t00, t01)/ max(t00, t01);
        final double altRatio = min(t10, t11)/ max(t10, t11);

        return Math.log(ratio) + Math.log(refRatio) - Math.log(altRatio);
    }

    /**
     * Returns an annotation result given a sor
     *
     * @param sor the symmetric odds ratio of the contingency table
     * @return a hash map from SOR
     */
    Map<String, Object> annotationForOneTable(final double sor) {
        return Collections.singletonMap(getKeyNames().get(0), formattedValue(sor));
    }

    public static String formattedValue(double sor) {
        return String.format("%.3f", sor);
    }


    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GaeaVCFConstants.STRAND_ODDS_RATIO_KEY);
    }
}
