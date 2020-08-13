package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public final class GenotypeSummaries extends InfoFieldAnnotation {

	@Override
	public List<String> getKeyNames() {
		return Arrays.asList(GaeaVCFConstants.NOCALL_CHROM_KEY, GaeaVCFConstants.GQ_MEAN_KEY,
				GaeaVCFConstants.GQ_STDEV_KEY);
	}

	@Override
	public Map<String, Object> annotate(ChromosomeInformationShare ref, VariantContext vc,
			ReadLikelihoods<Allele> likelihoods) {
		Utils.nonNull(vc);
		if (!vc.hasGenotypes()) {
			return Collections.emptyMap();
		}

		final Map<String, Object> returnMap = new LinkedHashMap<>();
		returnMap.put(GaeaVCFConstants.NOCALL_CHROM_KEY, vc.getNoCallCount());

		final DescriptiveStatistics stats = new DescriptiveStatistics();
		for (final Genotype g : vc.getGenotypes()) {
			if (g.hasGQ()) {
				stats.addValue(g.getGQ());
			}
		}
		if (stats.getN() > 0L) {
			returnMap.put(GaeaVCFConstants.GQ_MEAN_KEY, String.format("%.2f", stats.getMean()));
			if (stats.getN() > 1L) {
				returnMap.put(GaeaVCFConstants.GQ_STDEV_KEY, String.format("%.2f", stats.getStandardDeviation()));
			}
		}

		return returnMap;
	}
}
