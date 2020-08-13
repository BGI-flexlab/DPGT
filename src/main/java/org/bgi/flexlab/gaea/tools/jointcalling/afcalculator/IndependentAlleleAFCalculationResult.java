package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

import java.util.List;
import java.util.Map;

import htsjdk.variant.variantcontext.Allele;

class IndependentAlleleAFCalculationResult extends AFCalculationResult {
	/**
     * List of the supporting bi-allelic AFCalcResults that went into making this multi-allelic joint call
     */
    final List<AFCalculationResult> supporting;

    IndependentAlleleAFCalculationResult(final int[] alleleCountsOfMLE, final int nEvaluations,
                                         final List<Allele> allelesUsedInGenotyping, final double[] log10LikelihoodsOfAC,
                                         final double[] log10PriorsOfAC,
                                         final Map<Allele, Double> log10pRefByAllele, final List<AFCalculationResult> supporting) {
        super(alleleCountsOfMLE, nEvaluations, allelesUsedInGenotyping, log10LikelihoodsOfAC,
                log10PriorsOfAC, log10pRefByAllele);
        this.supporting = supporting;
    }
}
