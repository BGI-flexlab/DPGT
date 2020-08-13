package org.bgi.flexlab.gaea.tools.haplotypecaller.readfilter;

import java.io.Serializable;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.ReadUtils;

import htsjdk.samtools.SAMFileHeader;

public final class AlignmentAgreesWithHeaderReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1l;

    public AlignmentAgreesWithHeaderReadFilter( ) {}

    public AlignmentAgreesWithHeaderReadFilter( final SAMFileHeader header ) {
        super.setHeader(header);
    }

    @Override
    public boolean test( GaeaSamRecord read ) {
        return ReadUtils.alignmentAgreesWithHeader(samHeader, read);
    }
}
