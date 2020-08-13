package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.OptionalDouble;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.AlignmentUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

public final class ClippingRankSumTest extends RankSumTest implements StandardHCAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GaeaVCFConstants.CLIPPING_RANK_SUM_KEY); }

    @Override
    protected OptionalDouble getElementForRead(final GaeaSamRecord read, final int refLoc) {
        Utils.nonNull(read);
        return OptionalDouble.of(AlignmentUtils.getNumHardClippedBases(read));
    }
 }
