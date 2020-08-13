package org.bgi.flexlab.gaea.tools.haplotypecaller.smithwaterman;

public final class SWNativeAlignerResult {
    // CIGAR string of the alignment
    public final String cigar;

    // offset of the alignment
    public final int alignment_offset;


    public SWNativeAlignerResult(final String cigar, final int alignment_offset)
    {
        this.cigar = cigar;
        this.alignment_offset = alignment_offset;
    }
}
