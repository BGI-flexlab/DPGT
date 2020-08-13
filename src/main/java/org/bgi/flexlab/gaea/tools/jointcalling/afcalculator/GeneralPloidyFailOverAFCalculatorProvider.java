package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

public class GeneralPloidyFailOverAFCalculatorProvider extends AFCalculatorProvider {

    private final AFCalculator preferred;
    private final AFCalculatorImplementation preferredImplementation;
    private final AFCalculator failOver;

    /**
     * Creates a new AF calculator provider given the genotyping arguments and logger reference.
     * @param genotypeArgs genotyping parameter collection.
     * @param logger where the AF calculator logging messages go. If {@code null}, logging message will not be emitted.
     *
     * @throws IllegalStateException if {@code genotypeArgs} is {@code null}.
     */
    public GeneralPloidyFailOverAFCalculatorProvider(final int samplePloidy, int max_alternate_alleles) {
        preferredImplementation = AFCalculatorImplementation.bestValue(samplePloidy,max_alternate_alleles, null);
        preferred = preferredImplementation.newInstance();
        failOver = AFCalculatorImplementation.EXACT_GENERAL_INDEPENDENT.newInstance();
    }

    /**
     * {@inheritDoc}
     * @param ploidy {@inheritDoc}
     * @param maximumAlternativeAlleles {@inheritDoc}
     * @return {@inheritDoc}
     */
    @Override
    public AFCalculator getInstance(final int ploidy, final int maximumAlternativeAlleles) {
        return preferredImplementation.usableForParams(ploidy,maximumAlternativeAlleles) ? preferred : failOver;
    }
}

