package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

public class AS_BaseQualityRankSumTest extends AS_RankSumTest implements AS_StandardAnnotation {

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GaeaVCFConstants.AS_BASE_QUAL_RANK_SUM_KEY);
    }

    @Override
    public String getRawKeyName() { return GaeaVCFConstants.AS_RAW_BASE_QUAL_RANK_SUM_KEY;}

    /**
     * Get the element for the given read at the given reference position
     *
     * @param read     the read
     * @param refLoc   the reference position
     * @return a Double representing the element to be used in the rank sum test, or null if it should not be used
     */
    @Override
    protected OptionalDouble getElementForRead(final GaeaSamRecord read, final int refLoc) {
        return BaseQualityRankSumTest.getReadBaseQuality(read, refLoc);
    }

}

