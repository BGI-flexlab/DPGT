package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.haplotypecaller.annotator.HeterozygosityCalculator;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.AnnotationUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public final class AS_InbreedingCoeff extends InfoFieldAnnotation implements AS_StandardAnnotation {

    public static final int MIN_SAMPLES = 10;
    private Set<String> founderIds;    //TODO: either use this or enter a bug report

    public AS_InbreedingCoeff(){
        this(null);
    }

    public AS_InbreedingCoeff(final Set<String> founderIds){
        //If available, get the founder IDs and cache them. the IC will only be computed on founders then.
        this.founderIds = founderIds;
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GaeaVCFConstants.AS_INBREEDING_COEFFICIENT_KEY); }

    @VisibleForTesting
    public double calculateIC(final VariantContext vc, final Allele altAllele) {
        //make a new HeterozygosityUtils for each call to reset it
        return calculateIC(vc, altAllele, new HeterozygosityCalculator(vc));
    }

    private double calculateIC(final VariantContext vc, final Allele altAllele, final HeterozygosityCalculator heterozygosityUtils) {
        final int AN = vc.getCalledChrCount();
        final double altAF;

        final double hetCount = heterozygosityUtils.getHetCount(altAllele);

        final double F;
        //shortcut to get a value closer to the non-alleleSpecific value for bialleleics
        if (vc.isBiallelic()) {
            double refAC = heterozygosityUtils.getAlleleCount(vc.getReference());
            double altAC = heterozygosityUtils.getAlleleCount(altAllele);
            double refAF = refAC/(altAC+refAC);
            altAF = 1 - refAF;
            F = 1.0 - (hetCount / (2.0 * refAF * altAF * (double) heterozygosityUtils.getSampleCount())); // inbreeding coefficient
        }
        else {
            //compare number of hets for this allele (and any other second allele) with the expectation based on AFs
            //derive the altAF from the likelihoods to account for any accumulation of fractional counts from non-primary likelihoods,
            //e.g. for a GQ10 variant, the probability of the call will be ~0.9 and the second best call will be ~0.1 so adding up those 0.1s for het counts can dramatically change the AF compared with integer counts
            altAF = heterozygosityUtils.getAlleleCount(altAllele)/ (double) AN;
            F = 1.0 - (hetCount / (2.0 * (1 - altAF) * altAF * (double) heterozygosityUtils.getSampleCount())); // inbreeding coefficient
        }

        return F;
    }

	@Override
	public Map<String, Object> annotate(ChromosomeInformationShare ref, VariantContext vc,
			ReadLikelihoods<Allele> likelihoods) {
		Utils.nonNull(vc);
        final HeterozygosityCalculator heterozygosityUtils = new HeterozygosityCalculator(vc);

        if (heterozygosityUtils.getSampleCount() < MIN_SAMPLES) {
            return Collections.emptyMap();
        }
        final List<Double> ICvalues = new ArrayList<>();
        for (final Allele a : vc.getAlternateAlleles()) {
            ICvalues.add(calculateIC(vc, a, heterozygosityUtils));
        }
        return Collections.singletonMap(getKeyNames().get(0),  AnnotationUtils.encodeValueList(ICvalues, "%.4f"));
	}
}

