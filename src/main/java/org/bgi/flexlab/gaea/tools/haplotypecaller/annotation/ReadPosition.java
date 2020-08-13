package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.List;
import java.util.OptionalInt;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.AlignmentUtils;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.GaeaProtectedMathUtils;
import org.bgi.flexlab.gaea.util.ReadUtils;
import org.bgi.flexlab.gaea.util.Utils;

import com.google.common.primitives.Ints;

import htsjdk.variant.variantcontext.VariantContext;

public class ReadPosition extends PerAlleleAnnotation implements StandardMutectAnnotation {
    public static final String KEY = "MPOS";

    @Override
    protected int aggregate(final List<Integer> values) {
        return values.isEmpty() ? 0 : GaeaProtectedMathUtils.median(Ints.toArray(values));
    }

    @Override
    protected String getVcfKey() { return KEY; }

    @Override
    protected String getDescription() { return "median distance from end of read"; }

    @Override
    protected OptionalInt getValueForRead(final GaeaSamRecord read, final VariantContext vc) {
        Utils.nonNull(read);
        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), vc.getStart(), ReadUtils.ClippingTail.RIGHT_TAIL, true);
        if ( offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED || AlignmentUtils.isInsideDeletion(read.getCigar(), offset)) {
            return OptionalInt.empty();
        }

        int readPosition = AlignmentUtils.calcAlignmentByteArrayOffset(read.getCigar(), offset, false, 0, 0);
        final int numAlignedBases = AlignmentUtils.getNumAlignedBasesCountingSoftClips( read );
        final int distanceFromEnd = Math.min(readPosition, numAlignedBases - readPosition - 1);
        return OptionalInt.of(distanceFromEnd);
    }
}

