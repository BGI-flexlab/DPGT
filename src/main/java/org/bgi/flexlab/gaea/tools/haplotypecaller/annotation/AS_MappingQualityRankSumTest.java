package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

public class AS_MappingQualityRankSumTest extends AS_RankSumTest implements AS_StandardAnnotation {
    @Override
    public List<String> getKeyNames() { return Arrays.asList(GaeaVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY); }

    @Override
    public String getRawKeyName() { return GaeaVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY;}

    @Override
    protected OptionalDouble getElementForRead(final GaeaSamRecord read, final int refLoc) {
        Utils.nonNull(read);
        return OptionalDouble.of(read.getMappingQuality());
    }
}
