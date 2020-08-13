package org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation;

import java.util.List;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;

public abstract class GenotypeAnnotation extends VariantAnnotatorAnnotation {
	public abstract void annotate(final RefMetaDataTracker tracker, final ChromosomeInformationShare ref, final VariantContext vc,
			final Genotype g, final GenotypeBuilder gb);

	// return the descriptions used for the VCF FORMAT meta field
	public abstract List<VCFFormatHeaderLine> getDescriptions();
}
