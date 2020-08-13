package org.bgi.flexlab.gaea.tools.haplotypecaller.readfilter;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;

import htsjdk.samtools.SAMFileHeader;

public final class WellformedReadFilter extends ReadFilter {
    private static final long serialVersionUID = 1l;

    private ReadFilter wellFormedFilter = null;

    // Command line parser requires a no-arg constructor
    public WellformedReadFilter() {
    }

    @Override
    public void setHeader(SAMFileHeader header) {
        super.setHeader(header);
        createFilter();
    }

    public WellformedReadFilter(final SAMFileHeader header) {
        setHeader(header);
    }

    private void createFilter() {
        final AlignmentAgreesWithHeaderReadFilter alignmentAgreesWithHeader = new AlignmentAgreesWithHeaderReadFilter(samHeader);

        wellFormedFilter = ReadFilterLibrary.VALID_ALIGNMENT_START
                .and(ReadFilterLibrary.VALID_ALIGNMENT_END)
                .and(alignmentAgreesWithHeader)
                .and(ReadFilterLibrary.HAS_READ_GROUP)
                .and(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS)
                .and(ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH)
                .and(ReadFilterLibrary.SEQ_IS_STORED)
                .and(ReadFilterLibrary.CIGAR_CONTAINS_NO_N_OPERATOR);
    }

    @Override
    public boolean test(final GaeaSamRecord read ) {
        return wellFormedFilter.test(read);
    }
}
