package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public final class Coverage extends InfoFieldAnnotation implements StandardAnnotation, StandardMutectAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(VCFConstants.DEPTH_KEY); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

	@Override
	public Map<String, Object> annotate(ChromosomeInformationShare ref, VariantContext vc,
			ReadLikelihoods<Allele> likelihoods) {
		Utils.nonNull(vc);
        if (likelihoods == null || likelihoods.readCount() == 0) {
            return Collections.emptyMap();
        }

        final int depth = likelihoods.readCount();
        return Collections.singletonMap(getKeyNames().get(0), String.format("%d", depth));
	}
}

