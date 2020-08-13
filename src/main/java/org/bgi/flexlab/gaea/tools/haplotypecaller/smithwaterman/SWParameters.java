package org.bgi.flexlab.gaea.tools.haplotypecaller.smithwaterman;

public class SWParameters {
    private final int matchValue;
    private final int mismatchPenalty;
    private final int gapOpenPenalty;
    private final int gapExtendPenalty;

    /**
     * Create a new set of parameters for Smith-Waterman alignment
     * @param matchValue how much to reward a match during alignment >= 0
     * @param mismatchPenalty how much to penalize a mismatch during alignment <= 0
     * @param gapOpenPenalty how much penalize the creation of a new gap in the alignment <= 0
     * @param gapExtendPenalty how much to penalize extending an already open gap in the alignment <= 0
     */
    public SWParameters(final int matchValue, final int mismatchPenalty, final int gapOpenPenalty, final int gapExtendPenalty) {
        if( matchValue < 0 ) {
           throw new IllegalArgumentException("matchValue must be >= 0 but was passed as " + matchValue);
        }
        if( mismatchPenalty > 0 ) {
            throw new IllegalArgumentException("mismatchPenalty must be <= 0 but was passed as " + mismatchPenalty);
        }
        if( gapOpenPenalty > 0 ) {
            throw new IllegalArgumentException("gapOpenPenalty must be <= 0 but was passed as " + gapOpenPenalty);
        }
        if( gapExtendPenalty > 0 ) {
            throw new IllegalArgumentException("gapExtendPenalty must be <= 0 but was passed as " + gapExtendPenalty);
        }
        this.matchValue = matchValue;
        this.mismatchPenalty = mismatchPenalty;
        this.gapOpenPenalty = gapOpenPenalty;
        this.gapExtendPenalty = gapExtendPenalty;
    }

    /** gap extension penalty **/
    public int getGapExtendPenalty() {
        return gapExtendPenalty;
    }

    /** match value **/
    public int getMatchValue() {
        return matchValue;
    }

    /** mismatch penalty **/
    public int getMismatchPenalty() {
        return mismatchPenalty;
    }

    /** gap open penalty **/
    public int getGapOpenPenalty() {
        return gapOpenPenalty;
    }
}
