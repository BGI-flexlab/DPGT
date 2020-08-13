package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class HaplotypeScore extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation{

	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker, ChromosomeInformationShare ref, VariantContext vc) {
		return null;
	}

	@Override
    public List<String> getKeyNames() {
        return Arrays.asList(GaeaVCFConstants.HAPLOTYPE_SCORE_KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(GaeaVcfHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}
