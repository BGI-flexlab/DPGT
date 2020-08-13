package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.QualityUtils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public final class FisherStrand extends StrandBiasTest implements StandardAnnotation {

    static final double MIN_PVALUE = 1E-320;
    private static final int MIN_COUNT = ARRAY_DIM;

    // how large do we want the normalized table to be? (ie, sum of all entries must be smaller that this)
    private static final double TARGET_TABLE_SIZE = 200.0;

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GaeaVCFConstants.FISHER_STRAND_KEY);
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes){
        final int[][] tableFromPerSampleAnnotations = getTableFromSamples(genotypes, MIN_COUNT);
        return ( tableFromPerSampleAnnotations != null )? annotationForOneTable(pValueForContingencyTable(tableFromPerSampleAnnotations)) : null;
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromLikelihoods(final ReadLikelihoods<Allele> likelihoods,
                                                                     final VariantContext vc){
        final int[][] table = getContingencyTable(likelihoods, vc, MIN_COUNT);
        return annotationForOneTable(pValueForContingencyTable(table));
    }

    /**
     * Returns an annotation result given a pValue
     *
     * @param pValue
     * @return a hash map from FS -> phred-scaled pValue
     */
    Map<String, Object> annotationForOneTable(final double pValue) {
        return Collections.singletonMap(getKeyNames().get(0), makeValueObjectForAnnotation(pValue));
    }

    public static String makeValueObjectForAnnotation(final int[][] originalTable) {
        return makeValueObjectForAnnotation(pValueForContingencyTable(originalTable));
    }

    public static String makeValueObjectForAnnotation(double pValue) {
        return String.format("%.3f", QualityUtils.phredScaleErrorRate(Math.max(pValue, MIN_PVALUE))); // prevent INFINITYs
    }

    public static Double pValueForContingencyTable(final int[][] originalTable) {
        final int[][] normalizedTable = normalizeContingencyTable(originalTable);
        return FisherExactTest.twoSidedPValue(normalizedTable);
    }

    /**
     * Normalize the table so that the entries are not too large.
     * Note that this method does NOT necessarily make a copy of the table being passed in!
     *
     * @param table  the original table
     * @return a normalized version of the table or the original table if it is already normalized
     */
    private static int[][] normalizeContingencyTable(final int[][] table) {
        final int sum = addExact(table[0][0], table[0][1], table[1][0], table[1][1]);
        if ( sum <= TARGET_TABLE_SIZE * 2 ) {
            return table;
        }

        final double normFactor = sum / TARGET_TABLE_SIZE;

        return new int[][]{
                {(int) (table[0][0] / normFactor), (int) (table[0][1] / normFactor)},
                {(int) (table[1][0] / normFactor), (int) (table[1][1] / normFactor)},
        };
    }

    //Add a bunch of ints, blows up if there's overflow
    private static int addExact(final int... ints){
        int res = ints[0];
        for (int i = 1; i < ints.length; i++) {
            res = Math.addExact(res, ints[i]);
        }
        return res;
    }
}
