package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.GaeaVariantContextUtils;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public final class TandemRepeat extends InfoFieldAnnotation implements StandardMutectAnnotation {

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(
                GaeaVCFConstants.STR_PRESENT_KEY,
                GaeaVCFConstants.REPEAT_UNIT_KEY,
                GaeaVCFConstants.REPEATS_PER_ALLELE_KEY);
    }

    /*private static byte[] getRefBasesStartingAtVariantLocus(final ReferenceContext ref, final VariantContext vc) {
        final byte[] bases = ref.getBases();
        final int startIndex = vc.getStart() - ref.getWindow().getStart();
        return new String(bases).substring(startIndex).getBytes();
    }*/

	@Override
	public Map<String, Object> annotate(ChromosomeInformationShare ref, VariantContext vc,
			ReadLikelihoods<Allele> likelihoods) {
		/*Utils.nonNull(vc);
        if ( !vc.isIndel()) {
            return Collections.emptyMap();
        }

        final Pair<List<Integer>,byte[]> result = GaeaVariantContextUtils.getNumTandemRepeatUnits(vc, getRefBasesStartingAtVariantLocus(ref, vc));
        if (result == null) {
            return Collections.emptyMap();
        }

        final byte[] repeatUnit = result.getRight();
        final List<Integer> numUnits = result.getLeft();

        final Map<String, Object> map = new LinkedHashMap<>();
        map.put(GaeaVCFConstants.STR_PRESENT_KEY, true);
        map.put(GaeaVCFConstants.REPEAT_UNIT_KEY, new String(repeatUnit));
        map.put(GaeaVCFConstants.REPEATS_PER_ALLELE_KEY, numUnits);
        return Collections.unmodifiableMap(map);*/
		return null;
	}

}

