package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

import htsjdk.variant.variantcontext.VariantContext;

public final class FixedAFCalculatorProvider extends AFCalculatorProvider {

	private final AFCalculator calculator;

	private final boolean verifyRequests;

	private final int maximumAltAlleleCount;

	private final int ploidy;

	/**
	 * Constructs a fixed AF Calculator provider.
	 *
	 * @param configuration
	 *            the called configuration. This is the source of the fixed
	 *            ploidy and maximum number of supported alleles.
	 * @param verifyRequests
	 *            whether this provider will verify that each request for the AF
	 *            calculator meets the initial parameter values (ploidy,
	 *            sample-count and maximum number of alleles.
	 *
	 * @throws IllegalArgumentException
	 *             if {@code configuration} is {@code null}, or it contains
	 *             invalid values for sample ploidy and maximum number of
	 *             alternative alleles, or {@code sampleCount} is less than 0.
	 */
	public FixedAFCalculatorProvider(int samplePloidy, int max_alt_alleles, final boolean verifyRequests) {
		this(null, samplePloidy, max_alt_alleles, verifyRequests);
	}

	/**
	 * Constructs a fixed AF Calculator provider.
	 *
	 * @param preferred
	 *            preferred implementation.
	 * @param configuration
	 *            the called configuration. This is the source of the fixed
	 *            ploidy and maximum number of supported alleles.
	 * @param verifyRequests
	 *            whether this provider will verify that each request for the AF
	 *            calculator meets the initial parameter values (ploidy,
	 *            sample-count and maximum number of alleles.
	 *
	 * @throws IllegalArgumentException
	 *             if {@code configuration} is {@code null}, or it contains
	 *             invalid values for sample ploidy and maximum number of
	 *             alternative alleles, or {@code sampleCount} is less than 0.
	 */
	public FixedAFCalculatorProvider(final AFCalculatorImplementation preferred, int samplePloidy, int max_alt_alleles,
			final boolean verifyRequests) {
		if (samplePloidy < 1) {
			throw new IllegalArgumentException("invalid sample ploidy " + samplePloidy);
		}
		if (max_alt_alleles < 0) {
			throw new IllegalArgumentException("invalid maximum number of alleles " + (max_alt_alleles + 1));
		}

		ploidy = samplePloidy;
		maximumAltAlleleCount = max_alt_alleles;
		calculator = AFCalculatorImplementation.bestValue(ploidy, maximumAltAlleleCount, preferred).newInstance();
		this.verifyRequests = verifyRequests;
	}

	@Override
	public AFCalculator getInstance(final VariantContext vc, final int defaultPloidy, final int maximumAlleleCount) {
		if (verifyRequests) {
			// supers implementation will call eventually one of the other
			// methods, so no need to verify anything here.
			return super.getInstance(vc, defaultPloidy, maximumAlleleCount);
		}
		return calculator;
	}

	@Override
	public AFCalculator getInstance(final int ploidy, final int maxAltAlleleCount) {
		if (verifyRequests) {
			if (this.ploidy != AFCalculatorImplementation.UNBOUND_PLOIDY && ploidy != this.ploidy) {
				throw new IllegalArgumentException("non-supported ploidy:" + ploidy + " Only " + this.ploidy + " or "
						+ AFCalculatorImplementation.UNBOUND_PLOIDY);
			}
			if (maximumAltAlleleCount != AFCalculatorImplementation.UNBOUND_ALTERNATIVE_ALLELE_COUNT
					&& maxAltAlleleCount > maximumAltAlleleCount) {
				throw new IllegalArgumentException("non-supported alleleCount");
			}
		}
		return calculator;
	}

	/**
	 * Creates a fixed AF calculator provider that is thread safe.
	 *
	 * @param config
	 *            the caller configuration.
	 *
	 * @throws IllegalArgumentException
	 *             if any of the input argument is {@code null} or contain
	 *             invalid configuration like zero-samples, zero or negative
	 *             ploidy or negative-zero maximum number of alleles.
	 *
	 * @return never {@code null}
	 */
	public static AFCalculatorProvider createThreadSafeProvider(int samplePloidy, int max_alt_alleles) {

		return new ConcurrentAFCalculatorProvider() {
			@Override
			protected AFCalculatorProvider createProvider() {
				return new FixedAFCalculatorProvider(samplePloidy, max_alt_alleles, false);
			}
		};
	}
}
