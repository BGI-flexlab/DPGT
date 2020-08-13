package org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation;

import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.tools.jointcalling.annotator.AnnotationType;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public interface ActiveRegionBasedAnnotation extends AnnotationType {
    // return annotations for the given contexts split by sample and then read likelihood
    public abstract Map<String, Object> annotate(final VariantContext vc);

    // return the descriptions used for the VCF INFO meta field
    public abstract List<VCFInfoHeaderLine> getDescriptions();
}
