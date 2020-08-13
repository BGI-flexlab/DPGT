package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.bgi.flexlab.gaea.util.QualityUtils;

import cern.jet.math.Arithmetic;

public class StrandBiasTableUtils {

    //For now this is only for 2x2 contingency tables
    protected static final int ARRAY_DIM = 2;
    protected static final int ARRAY_SIZE = ARRAY_DIM * ARRAY_DIM;
    private static double MIN_PVALUE = 1E-320;
    // how large do we want the normalized table to be?
    private static final double TARGET_TABLE_SIZE = 200.0;
    private final static double AUGMENTATION_CONSTANT = 1.0;

    /**
     * Computes a two-sided p-Value for a Fisher's exact test on the contingency table, after normalizing counts so that the sum does not exceed {@value org.broadinstitute.gatk.tools.walkers.annotator.StrandBiasTableUtils#TARGET_TABLE_SIZE}
     * @param originalTable
     * @return
     */
    public static Double FisherExactPValueForContingencyTable(int[][] originalTable) {
    	final int[][] normalizedTable = normalizeContingencyTable(originalTable);

        int[][] table = copyContingencyTable(normalizedTable);

        int[] rowSums = { sumRow(table, 0), sumRow(table, 1) };
        int[] colSums = { sumColumn(table, 0), sumColumn(table, 1) };
        int N = rowSums[0] + rowSums[1];
        int sampleSize = colSums[0];
        int numberOfNonSuccesses = rowSums[1];
        int numberOfSuccesses = rowSums[0];
        /*
         * The lowest possible number of successes you can sample is what's left of your sample size after
         * drawing every non success in the urn. If the number of non successes in the urn is greater than the sample
         * size then the minimum number of drawn successes is 0.
         */
        int lo = Math.max(0, sampleSize - numberOfNonSuccesses);
        /*
         * The highest possible number of successes you can draw is either the total sample size or the number of
         * successes in the urn. (Whichever is smaller)
         */
        int hi = Math.min(sampleSize, numberOfSuccesses);

        /**
         * If the support of the distribution is only one value, creating the HypergeometricDistribution
         * doesn't make sense. There would be only one possible observation so the p-value has to be 1.
         */
        if (lo == hi) {
            return 1.0;
        }

        /**
         * For the hypergeometric distribution from which to calculate the probabilities:
         * The population (N = a+b+c+d) is the sum of all four numbers in the contingency table. Then the number of
         * "successes" (K = a+b) is the sum of the top row, and the sample size (n = a+c) is the sum of the first column.
         */
        final HypergeometricDistribution dist = new HypergeometricDistribution(N, numberOfSuccesses, sampleSize);

        //Then we determine a given probability with the sampled successes (k = a) from the first entry in the table.
        double pCutoff = dist.probability(table[0][0]);

        double pValue = 0.0;
        /**
         * Now cycle through all possible k's and add those whose probabilities are smaller than our observed k
         * to the p-value, since this is a two-sided test
         */
        for(int i = lo; i <= hi; i++){
            double pValuePiece = dist.probability(i);

            if(pValuePiece <= pCutoff) {
                pValue += pValuePiece;
            }
        }

        // min is necessary as numerical precision can result in pValue being slightly greater than 1.0
        return Math.min(pValue, 1.0);
    }

    /**
     * Helper function to turn the FisherStrand table into the SB annotation array
     * @param table the table used by the FisherStrand annotation
     * @return the array used by the per-sample Strand Bias annotation
     */
    public static List<Integer> getContingencyArray( final int[][] table ) {
        if(table.length != ARRAY_DIM || table[0].length != ARRAY_DIM) {
           return null;
        }

        final List<Integer> list = new ArrayList<>(ARRAY_SIZE);
        list.add(table[0][0]);
        list.add(table[0][1]);
        list.add(table[1][0]);
        list.add(table[1][1]);
        return list;
    }

    /**
     * Printing information to logger.info for debugging purposes
     *
     * @param name the name of the table
     * @param table the table itself
     */
    public static void printTable(final String name, final int[][] table) {
        final String pValue = String.format("%.3f", QualityUtils.phredScaleErrorRate(Math.max(FisherExactPValueForContingencyTable(table), MIN_PVALUE)));
    }

    /**
     * Adds the small value AUGMENTATION_CONSTANT to all the entries of the table.
     *
     * @param table the table to augment
     * @return the augmented table
     */
    public static double[][] augmentContingencyTable(final int[][] table) {
        double[][] augmentedTable = new double[ARRAY_DIM][ARRAY_DIM];
        for ( int i = 0; i < ARRAY_DIM; i++ ) {
            for ( int j = 0; j < ARRAY_DIM; j++ )
                augmentedTable[i][j] = table[i][j] + AUGMENTATION_CONSTANT;
        }

        return augmentedTable;
    }

    /**
     * Normalize the table so that the entries are not too large.
     * Note that this method does NOT necessarily make a copy of the table being passed in!
     *
     * @param table  the original table
     * @return a normalized version of the table or the original table if it is already normalized
     */
    public static int[][] normalizeContingencyTable(final int[][] table) {
        final int sum = table[0][0] + table[0][1] + table[1][0] + table[1][1];
        if ( sum <= TARGET_TABLE_SIZE * 2 )
            return table;

        final double normalizationFactor = (double)sum / TARGET_TABLE_SIZE;

        final int[][] normalized = new int[ARRAY_DIM][ARRAY_DIM];
        for ( int i = 0; i < ARRAY_DIM; i++ ) {
            for ( int j = 0; j < ARRAY_DIM; j++ )
                normalized[i][j] = (int)(table[i][j] / normalizationFactor);
        }

        return normalized;
    }

    public static int [][] copyContingencyTable(int [][] t) {
        int[][] c = new int[ARRAY_DIM][ARRAY_DIM];

        for ( int i = 0; i < ARRAY_DIM; i++ ) {
            //System.arraycopy(t,0,c,0,ARRAY_DIM);
            for (int j = 0; j < ARRAY_DIM; j++) {
                c[i][j] = t[i][j];
            }
        }

        return c;
    }

    public static boolean rotateTable(int[][] table) {
        table[0][0]--;
        table[1][0]++;

        table[0][1]++;
        table[1][1]--;

        return (table[0][0] >= 0 && table[1][1] >= 0);
    }

    public static boolean unrotateTable(int[][] table) {
        table[0][0]++;
        table[1][0]--;

        table[0][1]--;
        table[1][1]++;

        return (table[0][1] >= 0 && table[1][0] >= 0);
    }

    public static double computePValue(int[][] table) {

        int[] rowSums = { sumRow(table, 0), sumRow(table, 1) };
        int[] colSums = { sumColumn(table, 0), sumColumn(table, 1) };
        int N = rowSums[0] + rowSums[1];

        // calculate in log space for better precision
        double pCutoff = Arithmetic.logFactorial(rowSums[0])
                + Arithmetic.logFactorial(rowSums[1])
                + Arithmetic.logFactorial(colSums[0])
                + Arithmetic.logFactorial(colSums[1])
                - Arithmetic.logFactorial(table[0][0])
                - Arithmetic.logFactorial(table[0][1])
                - Arithmetic.logFactorial(table[1][0])
                - Arithmetic.logFactorial(table[1][1])
                - Arithmetic.logFactorial(N);
        return Math.exp(pCutoff);
    }

    public static int sumRow(int[][] table, int column) {
        int sum = 0;
        for (int r = 0; r < table.length; r++) {
            sum += table[r][column];
        }

        return sum;
    }

    private static int sumColumn(int[][] table, int row) {
        int sum = 0;
        for (int c = 0; c < table[row].length; c++) {
            sum += table[row][c];
        }

        return sum;
    }
}
