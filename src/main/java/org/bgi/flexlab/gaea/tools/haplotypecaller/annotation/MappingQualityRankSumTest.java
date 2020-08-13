package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.OptionalDouble;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

public final class MappingQualityRankSumTest extends RankSumTest implements StandardAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GaeaVCFConstants.MAP_QUAL_RANK_SUM_KEY); }

	@Override
	protected OptionalDouble getElementForRead(GaeaSamRecord read, int refLoc) {
		Utils.nonNull(read);
        return OptionalDouble.of(read.getMappingQuality());
	}
}
