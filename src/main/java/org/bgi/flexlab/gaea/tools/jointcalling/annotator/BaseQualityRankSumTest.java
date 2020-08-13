package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.Arrays;
import java.util.List;

import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class BaseQualityRankSumTest extends RankSumTest implements StandardAnnotation{
	@Override
    public List<String> getKeyNames() { return Arrays.asList(GaeaVCFConstants.BASE_QUAL_RANK_SUM_KEY); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(GaeaVcfHeaderLines.getInfoLine(getKeyNames().get(0))); }
}
