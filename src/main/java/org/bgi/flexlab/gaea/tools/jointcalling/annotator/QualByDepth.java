package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class QualByDepth extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation{
	
	protected final static double MAX_QD_BEFORE_FIXING = 35;
    protected final static double IDEAL_HIGH_QD = 30;
    protected final static double JITTER_SIGMA = 3;
    
	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker, ChromosomeInformationShare ref, VariantContext vc) {
		if ( !vc.hasLog10PError() )
            return null;

        final GenotypesContext genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        final int standardDepth = getDepth(genotypes);

        if ( standardDepth == 0 )
            return null;

        final double altAlleleLength = GaeaGvcfVariantContextUtils.getMeanAltAlleleLength(vc);
        	
        // Hack: UnifiedGenotyper (but not HaplotypeCaller or GenotypeGVCFs) over-estimates the quality of long indels
        //       Penalize the QD calculation for UG indels to compensate for this
        double QD = -10.0 * vc.getLog10PError() / ((double)standardDepth * indelNormalizationFactor(altAlleleLength, false));

        // Hack: see note in the fixTooHighQD method below
        QD = fixTooHighQD(QD);

        final Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), String.format("%.2f", QD));
        return map;
	}
	
	private double indelNormalizationFactor(final double altAlleleLength, final boolean increaseNormalizationAsLengthIncreases) {
        return ( increaseNormalizationAsLengthIncreases ? Math.max(altAlleleLength / 3.0, 1.0) : 1.0);
    }
	
	protected static double fixTooHighQD(final double QD) {
        if ( QD < MAX_QD_BEFORE_FIXING ) {
            return QD;
        } else {
        	double temp = GvcfMathUtils.getRandomGenerator().nextGaussian();
            return IDEAL_HIGH_QD + temp * JITTER_SIGMA;
        }
    }
	
	protected int getDepth(final GenotypesContext genotypes) {
        int standardDepth = 0;
        int ADrestrictedDepth = 0;

        for ( final Genotype genotype : genotypes ) {

            // we care only about variant calls with likelihoods
            if ( !genotype.isHet() && !genotype.isHomVar() )
                continue;

            // if we have the AD values for this sample, let's make sure that the variant depth is greater than 1!
            if ( genotype.hasAD() ) {
                final int[] AD = genotype.getAD();
                final int totalADdepth = (int)GvcfMathUtils.sum(AD);
                if ( totalADdepth - AD[0] > 1 )
                    ADrestrictedDepth += totalADdepth;
                standardDepth += totalADdepth;
                continue;
            }
            
            if ( genotype.hasDP() )
            	standardDepth += genotype.getDP();
        }

        // if the AD-restricted depth is a usable value (i.e. not zero), then we should use that one going forward
        if ( ADrestrictedDepth > 0 )
            standardDepth = ADrestrictedDepth;

        return standardDepth;
    }

	@Override
    public List<String> getKeyNames() { return Arrays.asList(GaeaVCFConstants.QUAL_BY_DEPTH_KEY); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(GaeaVcfHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

}
