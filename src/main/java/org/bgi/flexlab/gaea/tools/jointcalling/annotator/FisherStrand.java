package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;
import org.bgi.flexlab.gaea.tools.jointcalling.util.StrandBiasTableUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.QualityUtils;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class FisherStrand extends StrandBiasTest implements StandardAnnotation, ActiveRegionBasedAnnotation{

	private final static boolean ENABLE_DEBUGGING = false;

    private static final double MIN_PVALUE = 1E-320;
    private static final int MIN_QUAL_FOR_FILTERED_TEST = 17;
    private static final int MIN_COUNT = ARRAY_DIM;

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GaeaVCFConstants.FISHER_STRAND_KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(GaeaVcfHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes){
        final int[][] tableFromPerSampleAnnotations = getTableFromSamples( genotypes, MIN_COUNT );
        return ( tableFromPerSampleAnnotations != null )? pValueAnnotationForBestTable(tableFromPerSampleAnnotations, null) : null;
    }
    
    private Map<String, Object> pValueAnnotationForBestTable(final int[][] table1, final int[][] table2) {
        if ( table2 == null )
            return table1 == null ? null : annotationForOneTable(StrandBiasTableUtils.FisherExactPValueForContingencyTable(table1));
        else if (table1 == null)
            return annotationForOneTable(StrandBiasTableUtils.FisherExactPValueForContingencyTable(table2));
        else { // take the one with the best (i.e., least significant pvalue)
            double pvalue1 = StrandBiasTableUtils.FisherExactPValueForContingencyTable(table1);
            double pvalue2 = StrandBiasTableUtils.FisherExactPValueForContingencyTable(table2);
            return annotationForOneTable(Math.max(pvalue1, pvalue2));
        }
    }
    
    protected Map<String, Object> annotationForOneTable(final double pValue) {
        final Object value = String.format("%.3f", QualityUtils.phredScaleErrorRate(Math.max(pValue, MIN_PVALUE))); // prevent INFINITYs
        return Collections.singletonMap(getKeyNames().get(0), value);
    }
}
