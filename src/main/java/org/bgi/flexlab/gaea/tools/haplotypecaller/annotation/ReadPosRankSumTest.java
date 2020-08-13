package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.OptionalDouble;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.AlignmentUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.ReadUtils;
import org.bgi.flexlab.gaea.util.Utils;

public final class ReadPosRankSumTest extends RankSumTest implements StandardAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GaeaVCFConstants.READ_POS_RANK_SUM_KEY); }

    @Override
    protected OptionalDouble getElementForRead(final GaeaSamRecord read, final int refLoc) {
        return getReadPosition(read, refLoc);
    }

    @Override
    public boolean isUsableRead(final GaeaSamRecord read, final int refLoc) {
        Utils.nonNull(read);
        return super.isUsableRead(read, refLoc) && ReadUtils.getSoftEnd(read) >= refLoc;
    }

    public static OptionalDouble getReadPosition(final GaeaSamRecord read, final int refLoc) {
        Utils.nonNull(read);
        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), refLoc, ReadUtils.ClippingTail.RIGHT_TAIL, true);
        if ( offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
            return OptionalDouble.empty();
        }

        // If the offset inside a deletion, it does not lie on a read.
        if ( AlignmentUtils.isInsideDeletion(read.getCigar(), offset) ) {
            return OptionalDouble.of(INVALID_ELEMENT_FROM_READ);
        }

        int readPos = AlignmentUtils.calcAlignmentByteArrayOffset(read.getCigar(), offset, false, 0, 0);
        final int numAlignedBases = AlignmentUtils.getNumAlignedBasesCountingSoftClips( read );

        //After the middle of the read, we compute the postion from the end of the read.
        if (readPos > numAlignedBases / 2) {
            readPos = numAlignedBases - (readPos + 1);
        }
        return OptionalDouble.of(readPos);
    }


}
