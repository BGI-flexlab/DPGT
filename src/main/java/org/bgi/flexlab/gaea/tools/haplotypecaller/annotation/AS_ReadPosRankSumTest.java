package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.ReadUtils;
import org.bgi.flexlab.gaea.util.Utils;

public class AS_ReadPosRankSumTest extends AS_RankSumTest implements AS_StandardAnnotation {

    @Override
    public List<String> getKeyNames() { return Arrays.asList(GaeaVCFConstants.AS_READ_POS_RANK_SUM_KEY); }

    @Override
    public String getRawKeyName() { return GaeaVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY;}

    @Override
    protected OptionalDouble getElementForRead(final GaeaSamRecord read, final int refLoc) {
        return ReadPosRankSumTest.getReadPosition(read, refLoc);
    }

    @Override
    public boolean isUsableRead(final GaeaSamRecord read, final int refLoc) {
        Utils.nonNull(read);
        return super.isUsableRead(read, refLoc) && read.getSoftEnd() >= refLoc;
    }
}
