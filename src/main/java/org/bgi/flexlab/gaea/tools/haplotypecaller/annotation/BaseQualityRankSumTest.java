package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.OptionalDouble;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.ReadUtils;
import org.bgi.flexlab.gaea.util.Utils;

public final class BaseQualityRankSumTest extends RankSumTest implements StandardAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GaeaVCFConstants.BASE_QUAL_RANK_SUM_KEY); }

    @Override
    protected OptionalDouble getElementForRead(final GaeaSamRecord read, final int refLoc) {
        return getReadBaseQuality(read, refLoc);
    }

    public static OptionalDouble getReadBaseQuality(final GaeaSamRecord read, final int refLoc) {
        Utils.nonNull(read);
        return OptionalDouble.of(read.getBaseQuality(ReadUtils.getReadCoordinateForReferenceCoordinateUpToEndOfRead(read, refLoc, ReadUtils.ClippingTail.RIGHT_TAIL)));
    }
}
