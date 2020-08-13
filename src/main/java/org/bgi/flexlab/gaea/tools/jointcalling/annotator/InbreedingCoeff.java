package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.ReducibleAnnotationData;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ReducibleAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.AnnotationUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;
import org.bgi.flexlab.gaea.tools.jointcalling.util.HeterozygosityUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class InbreedingCoeff extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation, ReducibleAnnotation{

	protected static final int MIN_SAMPLES = 10;
    private Set<String> founderIds;
    private boolean didUniquifiedSampleNameCheck = false;
    protected HeterozygosityUtils heterozygosityUtils;
    final private boolean RETURN_ROUNDED = false;
    
    @Override
    public void initialize (Set<VCFHeaderLine> headerLines,Set<String> sampleList ) {
        //If available, get the founder IDs and cache them. the IC will only be computed on founders then.
        if(founderIds == null && sampleList != null) {
            founderIds = sampleList;
        }
        //intialize a HeterozygosityUtils before annotating for use in unit tests
        heterozygosityUtils = new HeterozygosityUtils(RETURN_ROUNDED);
    }
    
	@Override
	public String getRawKeyName() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Map<String, Object> annotateRawData(RefMetaDataTracker tracker, ChromosomeInformationShare ref,
			VariantContext vc) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Map<String, Object> combineRawData(List<Allele> allelesList,
			List<? extends ReducibleAnnotationData> listOfRawData) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Map<String, Object> finalizeRawData(VariantContext vc, VariantContext originalVC) {
		heterozygosityUtils = new HeterozygosityUtils(RETURN_ROUNDED);

        //if none of the "founders" are in the vc samples, assume we uniquified the samples upstream and they are all founders
        if (!didUniquifiedSampleNameCheck) {
            founderIds = AnnotationUtils.validateFounderIDs(founderIds, vc);
            didUniquifiedSampleNameCheck = true;
        }
        return makeCoeffAnnotation(vc);
	}

	@Override
	public void calculateRawData(VariantContext vc, ReducibleAnnotationData rawAnnotations) {
	}

	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker, ChromosomeInformationShare ref, VariantContext vc) {
		heterozygosityUtils = new HeterozygosityUtils(RETURN_ROUNDED);

        //if none of the "founders" are in the vc samples, assume we uniquified the samples upstream and they are all founders
        if (!didUniquifiedSampleNameCheck) {
            founderIds = AnnotationUtils.validateFounderIDs(founderIds, vc);
            didUniquifiedSampleNameCheck = true;
        }
        return makeCoeffAnnotation(vc);
	}
	
	protected Map<String, Object> makeCoeffAnnotation(final VariantContext vc) {
	    final GenotypesContext genotypes = (founderIds == null || founderIds.isEmpty()) ? vc.getGenotypes() : vc.getGenotypes(founderIds);
	    if (genotypes == null || genotypes.size() < MIN_SAMPLES || !vc.isVariant())
	        return null;
	    final double F = calculateIC(vc, genotypes);
	    if (heterozygosityUtils.getSampleCount() < MIN_SAMPLES)
	        return null;
	    return Collections.singletonMap(getKeyNames().get(0), (Object)String.format("%.4f", F));
	}
	
	protected double calculateIC(final VariantContext vc, final GenotypesContext genotypes) {

        final double[] genotypeCounts = heterozygosityUtils.getGenotypeCountsForRefVsAllAlts(vc, genotypes);  //guarantees that sampleCount is set
        if (genotypeCounts.length != 3) {
            throw new IllegalStateException("Input genotype counts must be length 3 for the number of genotypes with {2, 1, 0} ref alleles.");
        }
        final double refCount = genotypeCounts[HeterozygosityUtils.REF_INDEX];
        final double hetCount = genotypeCounts[HeterozygosityUtils.HET_INDEX];
        final double homCount = genotypeCounts[HeterozygosityUtils.VAR_INDEX];

        final double p = ( 2.0 * refCount + hetCount ) / ( 2.0 * (refCount + hetCount + homCount) ); // expected reference allele frequency
        final double q = 1.0 - p; // expected alternative allele frequency
        final double F = 1.0 - ( hetCount / ( 2.0 * p * q * (double) heterozygosityUtils.getSampleCount()) ); // inbreeding coefficient

        return F;
    }

	@Override
    public List<String> getKeyNames() { return Collections.singletonList(GaeaVCFConstants.INBREEDING_COEFFICIENT_KEY); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() { return Collections.singletonList(GaeaVcfHeaderLines.getInfoLine(getKeyNames().get(0))); }

}
