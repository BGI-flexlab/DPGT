package org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation;

import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.ReducibleAnnotationData;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.annotator.AnnotationType;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public interface ReducibleAnnotation extends AnnotationType{
	public abstract String getRawKeyName();

    /**
     * Generate the raw data necessary to calculate the annotation. Raw data is the final endpoint for gVCFs.
     *
     * @param tracker
     * @param ref
     * @param vc
     * @return
     */
    public abstract Map<String, Object> annotateRawData(final RefMetaDataTracker tracker,
                                                        final ChromosomeInformationShare ref,
                                                        final VariantContext vc);

    /**
     * Combine raw data, typically during the merging of raw data contained in multiple gVCFs as in CombineGVCFs and the
     * preliminary merge for GenotypeGVCFs
     * @param allelesList   The merged allele list across all variants being combined/merged
     * @param listOfRawData The raw data for all the variants being combined/merged
     * @return
     */
    public abstract Map<String, Object> combineRawData(final List<Allele> allelesList, final List <? extends ReducibleAnnotationData> listOfRawData);


    /**
     * Calculate the final annotation value from the raw data
     * @param vc -- contains the final set of alleles, possibly subset by GenotypeGVCFs
     * @param originalVC -- used to get all the alleles for all gVCFs
     * @return
     */
    public abstract Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC);

    /**
     *
     * @param vc
     * @param rawAnnotations
     */
    public abstract void calculateRawData(VariantContext vc, ReducibleAnnotationData rawAnnotations);
}
