package org.bgi.flexlab.gaea.tools.haplotypecaller.afcalculator;

import org.bgi.flexlab.gaea.tools.jointcalling.afcalculator.AFCalculator;
import org.bgi.flexlab.gaea.tools.jointcalling.afcalculator.AFCalculatorProvider;

import htsjdk.variant.variantcontext.VariantContext;

public abstract class ConcurrentAFCalculatorProvider extends AFCalculatorProvider {

    private final ThreadLocal<AFCalculatorProvider> threadLocal;

    /**
     * Create a new concurrent af-calculator provider instance.
     */
    public ConcurrentAFCalculatorProvider() {
        threadLocal = new ThreadLocal<AFCalculatorProvider>() {
            @Override
            public AFCalculatorProvider initialValue() {
                return createProvider();
            }
        };
    }

    @Override
    public AFCalculator getInstance(final VariantContext vc, final int defaultPloidy, final int maxAltAlleleCount) {
        return threadLocal.get().getInstance(vc,defaultPloidy,maxAltAlleleCount);
    }


    @Override
    public AFCalculator getInstance(final int ploidy, final int maxAltAlleleCount) {
        return threadLocal.get().getInstance(ploidy, maxAltAlleleCount);
    }

    protected abstract AFCalculatorProvider createProvider();
}

