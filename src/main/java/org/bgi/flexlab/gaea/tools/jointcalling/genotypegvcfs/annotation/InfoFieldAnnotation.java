package org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public abstract class InfoFieldAnnotation extends VariantAnnotatorAnnotation {

    public Map<String, Object> annotate(VariantContext vc) {
        return annotate(null, null, vc);
    }

    public Map<String, Object> annotate(ChromosomeInformationShare referenceContext, VariantContext vc) {

        return annotate(null, referenceContext, vc);
    }

    public abstract Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                                 final ChromosomeInformationShare ref,
                                                 final VariantContext vc);

    // return the descriptions used for the VCF INFO meta field
    public List<VCFInfoHeaderLine> getDescriptions() {
        final List<VCFInfoHeaderLine> lines = new ArrayList<>(5);
        for (final String key : getKeyNames()) {
            lines.add(GaeaVcfHeaderLines.getInfoLine(key));
        }
        return lines;
    }
}
