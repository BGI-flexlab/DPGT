package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ReducibleAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class RMSMappingQuality extends RMSAnnotation
		implements StandardAnnotation, ActiveRegionBasedAnnotation, ReducibleAnnotation {

	@Override
	public String getRawKeyName() {
		return GaeaVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY;
	}

	@Override
	protected String makeRawAnnotationString(List<Allele> vcAlleles, Map<Allele, Number> sumOfSquares) {
		return String.format("%.2f", sumOfSquares.get(Allele.NO_CALL));
	}

	@Override
	protected String makeFinalizedAnnotationString(VariantContext vc, Map<Allele, Number> sumOfSquares) {
		int numOfReads = getNumOfReads(vc);
        return String.format("%.2f", Math.sqrt((double)sumOfSquares.get(Allele.NO_CALL)/numOfReads));
	}

	@Override
	public List<VCFInfoHeaderLine> getDescriptions() {
		final List<VCFInfoHeaderLine> headerLines = new ArrayList<>();
		headerLines.add(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
		return headerLines;
	}

	public List<String> getKeyNames() {
		return Arrays.asList(VCFConstants.RMS_MAPPING_QUALITY_KEY);
	}
	
	
}
