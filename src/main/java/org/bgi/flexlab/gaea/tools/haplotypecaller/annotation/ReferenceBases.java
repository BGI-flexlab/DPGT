package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class ReferenceBases extends InfoFieldAnnotation {
    public static final String REFERENCE_BASES_KEY = "REF_BASES";

    public static final int NUM_BASES_ON_EITHER_SIDE = 10;

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(REFERENCE_BASES_KEY); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(ReferenceBases.REFERENCE_BASES_KEY, 1, VCFHeaderLineType.String, "local reference bases."));
    }

	@Override
	public Map<String, Object> annotate(ChromosomeInformationShare ref, VariantContext vc,
			ReadLikelihoods<Allele> likelihoods) {
		return null;
		/*if (ref==null)  {
            return Collections.emptyMap();
        }
        final int basesToDiscardInFront = Math.max(vc.getStart() - ref.getWindow().getStart() - NUM_BASES_ON_EITHER_SIDE, 0);
        final String allBases = new String(ref.getBases());
        final String localBases = allBases.substring(basesToDiscardInFront, basesToDiscardInFront + 2 * NUM_BASES_ON_EITHER_SIDE);
        return Collections.singletonMap(REFERENCE_BASES_KEY, localBases );*/
	}
}

