package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public final class SampleList extends InfoFieldAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GaeaVCFConstants.SAMPLE_LIST_KEY); }

	@Override
	public Map<String, Object> annotate(ChromosomeInformationShare ref, VariantContext vc,
			ReadLikelihoods<Allele> likelihoods) {
		Utils.nonNull(vc);
        if ( vc.isMonomorphicInSamples() || !vc.hasGenotypes() ) {
            return Collections.emptyMap();
        }

        final StringBuilder samples = new StringBuilder();
        for ( final Genotype genotype : vc.getGenotypesOrderedByName() ) {
            if ( genotype.isCalled() && !genotype.isHomRef() ){
                if ( samples.length() > 0 ) {
                    samples.append(",");
                }
                samples.append(genotype.getSampleName());
            }
        }

        if ( samples.length() == 0 ) {
            return Collections.emptyMap();
        }

        return Collections.singletonMap(getKeyNames().get(0), samples.toString());
	}
}

