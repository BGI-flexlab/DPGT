package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.Utils;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public final class MappingQualityZero extends InfoFieldAnnotation {

    @VisibleForTesting
    static String formattedValue(long mq0) {
        return String.format("%d", mq0);
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(VCFConstants.MAPPING_QUALITY_ZERO_KEY); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

	@Override
	public Map<String, Object> annotate(ChromosomeInformationShare ref, VariantContext vc,
			ReadLikelihoods<Allele> likelihoods) {
		Utils.nonNull(vc);
        if (!vc.isVariant() || likelihoods == null){
            return Collections.emptyMap();
        }
        //NOTE: unlike other annotations, this one returns 0 if likelihoods are empty
        final long mq0 = IntStream.range(0, likelihoods.numberOfSamples()).boxed()
                .flatMap(s -> likelihoods.sampleReads(s).stream())
                .filter(r -> r.getMappingQuality() == 0)
                .count();

        return Collections.singletonMap(getKeyNames().get(0), formattedValue(mq0));
	}
}

