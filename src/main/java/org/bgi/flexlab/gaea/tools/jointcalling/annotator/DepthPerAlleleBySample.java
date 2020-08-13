package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.Arrays;
import java.util.List;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.GenotypeAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class DepthPerAlleleBySample extends GenotypeAnnotation implements StandardAnnotation{

	@Override
	public void annotate(RefMetaDataTracker tracker, ChromosomeInformationShare ref, VariantContext vc, Genotype g,
			GenotypeBuilder gb) {
		if ( g == null || !g.isCalled() )
            return;
	}

	@Override
	public List<VCFFormatHeaderLine> getDescriptions() {
		return Arrays.asList(VCFStandardHeaderLines.getFormatLine(getKeyNames().get(0)));
	}

	@Override
	public List<String> getKeyNames() {
		return Arrays.asList(VCFConstants.GENOTYPE_ALLELE_DEPTHS);
	}

}
