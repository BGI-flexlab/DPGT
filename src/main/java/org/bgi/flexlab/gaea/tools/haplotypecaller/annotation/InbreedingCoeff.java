package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.GenotypeCounts;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.GenotypeUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public final class InbreedingCoeff extends InfoFieldAnnotation implements StandardAnnotation {

    private static final Logger logger = LogManager.getLogger(InbreedingCoeff.class);
    private static final int MIN_SAMPLES = 10;
    private static final boolean ROUND_GENOTYPE_COUNTS = false;
    private final Set<String> founderIds;

    public InbreedingCoeff(){
        this(null);
    }

    public InbreedingCoeff(final Set<String> founderIds){
        //If available, get the founder IDs and cache them. the IC will only be computed on founders then.
        this.founderIds = founderIds;
    }

    @VisibleForTesting
    static Pair<Integer, Double> calculateIC(final VariantContext vc, final GenotypesContext genotypes) {
        final GenotypeCounts t = GenotypeUtils.computeDiploidGenotypeCounts(vc, genotypes, ROUND_GENOTYPE_COUNTS);

        final double refCount = t.getRefs();
        final double hetCount = t.getHets();
        final double homCount = t.getHoms();
        // number of samples that have likelihoods
        final int sampleCount = (int) genotypes.stream().filter(g-> GenotypeUtils.isDiploidWithLikelihoods(g)).count();

        final double p = ( 2.0 * refCount + hetCount ) / ( 2.0 * (refCount + hetCount + homCount) ); // expected reference allele frequency
        final double q = 1.0 - p; // expected alternative allele frequency
        final double expectedHets = 2.0 * p * q * sampleCount; //numbers of hets that would be expected based on the allele frequency (asuming Hardy Weinberg Equilibrium)
        final double F = 1.0 - ( hetCount / expectedHets ); // inbreeding coefficient

        return Pair.of(sampleCount, F);
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GaeaVCFConstants.INBREEDING_COEFFICIENT_KEY); }

	@Override
	public Map<String, Object> annotate(ChromosomeInformationShare ref, VariantContext vc,
			ReadLikelihoods<Allele> likelihoods) {
		Utils.nonNull(vc);
        final GenotypesContext genotypes = (founderIds == null || founderIds.isEmpty()) ? vc.getGenotypes() : vc.getGenotypes(founderIds);
        if (genotypes == null || genotypes.size() < MIN_SAMPLES || !vc.isVariant()) {
            return Collections.emptyMap();
        }
        final Pair<Integer, Double> sampleCountCoeff = calculateIC(vc, genotypes);
        final int sampleCount = sampleCountCoeff.getLeft();
        final double F = sampleCountCoeff.getRight();
        if (sampleCount < MIN_SAMPLES) {
            logger.warn("Annotation will not be calculated, must provide at least " + MIN_SAMPLES + " samples");
            return Collections.emptyMap();
        }
        return Collections.singletonMap(getKeyNames().get(0), String.format("%.4f", F));
	}

}

