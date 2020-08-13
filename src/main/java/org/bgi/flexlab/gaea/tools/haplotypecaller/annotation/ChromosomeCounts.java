package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public final class ChromosomeCounts extends InfoFieldAnnotation implements StandardAnnotation {

    public static final String[] keyNames = {
            VCFConstants.ALLELE_NUMBER_KEY,
            VCFConstants.ALLELE_COUNT_KEY,
            VCFConstants.ALLELE_FREQUENCY_KEY };

    public static final VCFInfoHeaderLine[] descriptions = {
            VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY),
            VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY),
            VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY) };

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(keyNames);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(descriptions); }

	@Override
	public Map<String, Object> annotate(ChromosomeInformationShare ref, VariantContext vc,
			ReadLikelihoods<Allele> likelihoods) {
		Utils.nonNull(vc);
        if ( ! vc.hasGenotypes() ) {
            return Collections.emptyMap();
        }

        return VariantContextUtils.calculateChromosomeCounts(vc, new LinkedHashMap<>(), true, Collections.emptySet());
	}
}

