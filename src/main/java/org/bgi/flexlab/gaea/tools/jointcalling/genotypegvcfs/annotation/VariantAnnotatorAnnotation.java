package org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation;

import java.util.List;
import java.util.Set;

import htsjdk.variant.vcf.VCFHeaderLine;

public abstract class VariantAnnotatorAnnotation {
	// return the INFO keys
    public abstract List<String> getKeyNames();

    // initialization method (optional for subclasses, and therefore non-abstract)
    public void initialize (Set<VCFHeaderLine> headerLines,Set<String> sampleList) { }
}
