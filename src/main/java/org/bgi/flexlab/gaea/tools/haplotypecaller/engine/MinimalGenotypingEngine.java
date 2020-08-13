package org.bgi.flexlab.gaea.tools.haplotypecaller.engine;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.haplotypecaller.SampleList;
import org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection.UnifiedArgumentCollection;
import org.bgi.flexlab.gaea.tools.jointcalling.UnifiedGenotypingEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.UnifiedGenotypingEngine.GenotypingOutputMode;
import org.bgi.flexlab.gaea.tools.jointcalling.UnifiedGenotypingEngine.OutputMode;
import org.bgi.flexlab.gaea.tools.jointcalling.afcalculator.AFCalculatorProvider;

import htsjdk.variant.variantcontext.Allele;

public final class MinimalGenotypingEngine extends GenotypingEngine<UnifiedArgumentCollection> {

    /**
     * Creates a new unified genotyping given the UG configuration parameters and the targeted set of samples
     *
     * @param configuration the UG configuration.
     * @param samples list of samples
     */
    public MinimalGenotypingEngine(final UnifiedArgumentCollection configuration, final SampleList samples,
                                   final AFCalculatorProvider afCalculatorProvider) {
        this(configuration, samples, afCalculatorProvider, false);
    }

    /**
     * Creates a new unified genotyping given the UG configuration parameters and the targeted set of samples
     *
     * @param configuration the UG configuration.
     * @param samples list of samples
     * @param doAlleleSpecificCalcs Whether to calculate genotyping annotations needed for allele specific annotations
     */
    public MinimalGenotypingEngine(final UnifiedArgumentCollection configuration, final SampleList samples,
                                    final AFCalculatorProvider afCalculatorProvider, boolean doAlleleSpecificCalcs ) {
        super(configuration, samples, afCalculatorProvider, doAlleleSpecificCalcs);

        if ( configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            throw new UserException("GENOTYPE_GIVEN_ALLELES mode not supported in the MinimalGenotypingEngine");
        } else if ( configuration.GLmodel != UnifiedGenotypingEngine.Model.SNP ) {
            throw new UserException("Only the diploid SNP model is supported in the MinimalGenotypingEngine");
        } else if ( configuration.COMPUTE_SLOD ) {
            throw new UserException("--computeSLOD not supported in the MinimalGenotypingEngine");
        }
    }

    @Override
    protected boolean forceKeepAllele(final Allele allele) {
        return configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES || configuration.annotateAllSitesWithPLs;
    }

    @Override
    protected String callSourceString() {
        return "UG_call";
    }

    @Override
    protected boolean forceSiteEmission() {
        return configuration.outputMode == OutputMode.EMIT_ALL_SITES;
    }
}
