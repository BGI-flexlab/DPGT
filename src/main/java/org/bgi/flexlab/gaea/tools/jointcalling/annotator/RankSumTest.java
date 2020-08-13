package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public class RankSumTest extends InfoFieldAnnotation implements ActiveRegionBasedAnnotation{
	protected static double INVALID_ELEMENT_FROM_READ = Double.NEGATIVE_INFINITY;

	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker, ChromosomeInformationShare ref, VariantContext vc) {
		final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || genotypes.isEmpty())
            return null;

        final List<Double> refQuals = new ArrayList<>();
        final List<Double> altQuals = new ArrayList<>();

        if ( refQuals.isEmpty() && altQuals.isEmpty() )
            return null;

        final MannWhitneyU mannWhitneyU = new MannWhitneyU();
        
        // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
        final MannWhitneyU.Result result = mannWhitneyU.test(convertToArray(altQuals), convertToArray(refQuals), MannWhitneyU.TestType.FIRST_DOMINATES);
        final double zScore = result.getZ();


        final Map<String, Object> map = new HashMap<>();
        if (!Double.isNaN(zScore))
            map.put(getKeyNames().get(0), String.format("%.3f", zScore));
        return map;
	}
	
	public static double[] convertToArray(List<Double> list){
        double[] ret = new double[list.size()];
        Iterator<Double> iterator = list.iterator();
        for (int i = 0; i < ret.length; i++)
        {
            ret[i] = iterator.next().doubleValue();
        }
        return ret;
    }

	@Override
	public List<String> getKeyNames() {
		return null;
	}

}
