package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.ChromosomeCountConstants;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class ChromosomeCounts extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation{

	private Set<String> founderIds = new HashSet<String>();
	
	private boolean didUniquifiedSampleNameCheck = false;
	
	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker, ChromosomeInformationShare ref, VariantContext vc) {
		if ( ! vc.hasGenotypes() )
            return null;
        //if none of the "founders" are in the vc samples, assume we uniquified the samples upstream and they are all founders
        if (!didUniquifiedSampleNameCheck) {
            checkSampleNames(vc);
            didUniquifiedSampleNameCheck = true;
        }

        return VariantContextUtils.calculateChromosomeCounts(vc, new HashMap<String, Object>(), true,founderIds);
	}
	
	@Override
	public void initialize ( Set<VCFHeaderLine> headelines,Set<String> sampleList ){
        //If families were given, get the founders ids
        founderIds = sampleList;
    }

	public List<String> getKeyNames() {
        return Arrays.asList(ChromosomeCountConstants.keyNames);
    }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(ChromosomeCountConstants.descriptions); }

    //this method is intended to reconcile uniquified sample names
    // it comes into play when calling this annotation from GenotypeGVCFs with --uniquifySamples because founderIds
    // is derived from the sampleDB, which comes from the input sample names, but vc will have uniquified (i.e. different)
    // sample names. Without this check, the founderIds won't be found in the vc and the annotation won't be calculated.
    protected void checkSampleNames(final VariantContext vc) {
        Set<String> vcSamples = new HashSet<>();
        vcSamples.addAll(vc.getSampleNames());
        if (!vcSamples.isEmpty()) {
            if (founderIds!=null) {
                vcSamples.retainAll(founderIds);
                if (vcSamples.isEmpty())
                    founderIds = vc.getSampleNames();
            }
        }
    }

}
