package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;
import org.bgi.flexlab.gaea.tools.jointcalling.util.StrandBiasTableUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class StrandOddsRatio extends StrandBiasTest implements StandardAnnotation, ActiveRegionBasedAnnotation{
	private static final int MIN_COUNT = 0;
	
	final protected double calculateSOR(final int[][] originalTable) {
	    final double[][] augmentedTable = StrandBiasTableUtils.augmentContingencyTable(originalTable);

	    double ratio = 0;

	    ratio += (augmentedTable[0][0] / augmentedTable[0][1]) * (augmentedTable[1][1] / augmentedTable[1][0]);
	    ratio += (augmentedTable[0][1] / augmentedTable[0][0]) * (augmentedTable[1][0] / augmentedTable[1][1]);

	    final double refRatio = (Math.min(augmentedTable[0][0], augmentedTable[0][1])/Math.max(augmentedTable[0][0], augmentedTable[0][1]));
	    final double altRatio = (Math.min(augmentedTable[1][0], augmentedTable[1][1])/Math.max(augmentedTable[1][0], augmentedTable[1][1]));

	    ratio = ratio*refRatio/altRatio;

	    return Math.log(ratio);
	}
	
	protected Map<String, Object> annotationForOneTable(final double ratio) {
	    final Object value = String.format("%.3f", ratio);
	    return Collections.singletonMap(getKeyNames().get(0), value);
	}

	@Override
	protected Map<String, Object> calculateAnnotationFromGTfield(GenotypesContext genotypes) {
		final int[][] tableFromPerSampleAnnotations = getTableFromSamples( genotypes, MIN_COUNT );
        if ( tableFromPerSampleAnnotations != null ) {
            final double ratio = calculateSOR(tableFromPerSampleAnnotations);
            return annotationForOneTable(ratio);
        }
        return null;
	}

	@Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(GaeaVcfHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GaeaVCFConstants.STRAND_ODDS_RATIO_KEY);
    }

}
