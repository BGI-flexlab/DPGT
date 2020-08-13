package org.bgi.flexlab.gaea.tools.haplotypecaller.utils;

import java.io.Serializable;
import java.util.Comparator;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.ReadUtils;

import htsjdk.samtools.SAMFileHeader;

public final class ReadCoordinateComparator implements Comparator<GaeaSamRecord>, Serializable {
    private static final long serialVersionUID = 1L;

    private final SAMFileHeader header;

    public ReadCoordinateComparator( final SAMFileHeader header ) {
        if ( header == null ) {
            throw new IllegalArgumentException("header must be non-null");
        }
        this.header = header;
    }

    @Override
    public int compare( GaeaSamRecord first, GaeaSamRecord second ) {
        int result = compareCoordinates(first, second, header);
        if ( result != 0 ) {
            return result;
        }

        //This is done to mimic SAMRecordCoordinateComparator's behavior
        if (first.getReadNegativeStrandFlag() != second.getReadNegativeStrandFlag()) {
            return first.getReadNegativeStrandFlag()? 1: -1;
        }

        if ( first.getReadName() != null && second.getReadName() != null ) {
            result = first.getReadName().compareTo(second.getReadName());
            if ( result != 0 ) { return result; }
        }
        result = Integer.compare(ReadUtils.getSAMFlagsForRead(first), ReadUtils.getSAMFlagsForRead(second));
        if ( result != 0 ) { return result; }
        result = Integer.compare(first.getMappingQuality(), second.getMappingQuality());
        if ( result != 0 ) { return result; }
        if (first.getReadPairedFlag() && second.getReadPairedFlag()) {
            result = Integer.compare(ReadUtils.getMateReferenceIndex(first, header), ReadUtils.getMateReferenceIndex(second, header));
            if ( result != 0 ) { return result; }
            result = Integer.compare(first.getMateStart(), second.getMateStart());
            if ( result != 0 ) { return result; }
        }
        result = Integer.compare(first.getInferredInsertSize(), second.getInferredInsertSize());

        return result;
    }

    public static int compareCoordinates( final GaeaSamRecord first, final GaeaSamRecord second, final SAMFileHeader header ) {
        final int firstRefIndex = ReadUtils.getAssignedReferenceIndex(first, header);
        final int secondRefIndex = ReadUtils.getAssignedReferenceIndex(second, header);

        if ( firstRefIndex == -1 ) {
            return (secondRefIndex == -1 ? 0 : 1);
        }
        else if ( secondRefIndex == -1 ) {
            return -1;
        }

        final int refIndexDifference = firstRefIndex - secondRefIndex;
        if ( refIndexDifference != 0 ) {
            return refIndexDifference;
        }

        return Integer.compare(first.getAlignmentStart(), second.getAlignmentStart());
    }
}
