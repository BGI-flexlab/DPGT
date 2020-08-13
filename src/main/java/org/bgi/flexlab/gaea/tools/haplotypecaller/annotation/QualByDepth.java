package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public final class QualByDepth extends InfoFieldAnnotation implements StandardAnnotation {

    static final double MAX_QD_BEFORE_FIXING = 35;

    static final double IDEAL_HIGH_QD = 30;
    private static final double JITTER_SIGMA = 3;

    @Override
    public Map<String, Object> annotate(final ChromosomeInformationShare ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if ( !vc.hasLog10PError() ) {
            return Collections.emptyMap();
        }

        final GenotypesContext genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.isEmpty() ) {
            return Collections.emptyMap();
        }

        int depth = 0;
        int ADrestrictedDepth = 0;

        for ( final Genotype genotype : genotypes ) {
            // we care only about variant calls with likelihoods
            if ( !genotype.isHet() && !genotype.isHomVar() ) {
                continue;
            }

            // if we have the AD values for this sample, let's make sure that the variant depth is greater than 1!
            if ( genotype.hasAD() ) {
                final int[] AD = genotype.getAD();
                final int totalADdepth = (int) GvcfMathUtils.sum(AD);
                if ( totalADdepth != 0 ) {
                    if (totalADdepth - AD[0] > 1) {
                        ADrestrictedDepth += totalADdepth;
                    }
                    depth += totalADdepth;
                }
            } else if (likelihoods != null) {
                depth += likelihoods.sampleReadCount(likelihoods.indexOfSample(genotype.getSampleName()));
            } else if ( genotype.hasDP() ) {
                depth += genotype.getDP();
            }
        }

        // if the AD-restricted depth is a usable value (i.e. not zero), then we should use that one going forward
        if ( ADrestrictedDepth > 0 ) {
            depth = ADrestrictedDepth;
        }

        if ( depth == 0 ) {
            return Collections.emptyMap();
        }

        final double qual = -10.0 * vc.getLog10PError();
        double QD = qual / depth;

        // Hack: see note in the fixTooHighQD method below
        QD = fixTooHighQD(QD);

        return Collections.singletonMap(getKeyNames().get(0), String.format("%.2f", QD));
    }

    /**
     * The haplotype caller generates very high quality scores when multiple events are on the
     * same haplotype.  This causes some very good variants to have unusually high QD values,
     * and VQSR will filter these out.  This code looks at the QD value, and if it is above
     * threshold we map it down to the mean high QD value, with some jittering
     *
     * @param QD the raw QD score
     * @return a QD value
     */
    public static double fixTooHighQD(final double QD) {
        if ( QD < MAX_QD_BEFORE_FIXING ) {
            return QD;
        } else {
            return IDEAL_HIGH_QD + GvcfMathUtils.getRandomGenerator().nextGaussian() * JITTER_SIGMA;
        }
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GaeaVCFConstants.QUAL_BY_DEPTH_KEY); }
}
