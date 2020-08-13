package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

import java.lang.reflect.Constructor;
import java.lang.reflect.Modifier;
import java.util.HashMap;
import java.util.Map;

public enum AFCalculatorImplementation {
	/**
	 * Fast implementation for multi-allelics (equivalent to
	 * {@link #EXACT_REFERENCE} for biallelics sites
	 */
	EXACT_INDEPENDENT(IndependentAllelesDiploidExactAFCalculator.class, 2),

	/**
	 * reference implementation of multi-allelic EXACT model. Extremely slow for
	 * many alternate alleles
	 */
	EXACT_REFERENCE(ReferenceDiploidExactAFCalculator.class, 2),

	/** original biallelic exact model, for testing only */
	EXACT_ORIGINAL(OriginalDiploidExactAFCalculator.class, 2, 2),

	/**
	 * implementation that supports any sample ploidy. Currently not available
	 * for the HaplotypeCaller
	 */
	EXACT_GENERAL_PLOIDY(GeneralPloidyExactAFCalculator.class),

	/**
	 * Implementation that implements the {@link #EXACT_INDEPENDENT} for any
	 * ploidy.
	 */
	EXACT_GENERAL_INDEPENDENT(IndependentAllelesExactAFCalculator.class);

	/**
	 * Special max alt allele count indicating that this maximum is in fact
	 * unbound (can be anything).
	 */
	public final static int UNBOUND_ALTERNATIVE_ALLELE_COUNT = -1;

	/**
	 * Special ploidy constant that indicates that in fact the ploidy is unbound
	 * (can be anything).
	 */
	public final static int UNBOUND_PLOIDY = -1;

	private static Map<Class<? extends AFCalculator>, AFCalculatorImplementation> calculatorClassToValue = buildCalculatorClassToValueMap();

	/**
	 * Reference to the calculator class.
	 */
	public final Class<? extends AFCalculator> calculatorClass;

	/**
	 * Maximum number of supported alternative alleles.
	 */
	public final int maxAltAlleles;

	/**
	 * Reference to the constructor to instantiate a calculator for this
	 * implementation.
	 */
	protected final Constructor<? extends AFCalculator> constructor;

	/**
	 * Supported ploidy.
	 *
	 * This is equal to {@link #UNBOUND_PLOIDY} if the class can handle any
	 * ploidy.
	 */
	public final int requiredPloidy;

	/**
	 * Reference to the default implementation.
	 */
	public final static AFCalculatorImplementation DEFAULT = EXACT_INDEPENDENT;

	AFCalculatorImplementation(final Class<? extends AFCalculator> clazz, final int requiredPloidy,
			final int maxAltAlleles) {
		calculatorClass = clazz;
		this.requiredPloidy = requiredPloidy;
		this.maxAltAlleles = maxAltAlleles;
		this.constructor = findInstantiationConstructor(calculatorClass);
	}

	AFCalculatorImplementation(final Class<? extends AFCalculator> clazz) {
		this(clazz, UNBOUND_PLOIDY, UNBOUND_ALTERNATIVE_ALLELE_COUNT);
	}

	/**
	 * Constructs a new instance leaving max-allele count unbound.
	 * 
	 * @param clazz
	 *            the calculator class that realizes this implementation.
	 * @param requiredPloidy
	 *            the required ploidy; zero or greater or
	 *            {@link #UNBOUND_PLOIDY} to indicate that any ploidy is
	 *            supported.
	 */
	AFCalculatorImplementation(final Class<? extends AFCalculator> clazz, final int requiredPloidy) {
		this(clazz, requiredPloidy, UNBOUND_PLOIDY);
	}

	/**
	 * Checks whether a given ploidy and max alternative alleles combination is
	 * supported or not.
	 * 
	 * @param requestedPloidy
	 *            the targeted ploidy.
	 * @param requestedMaxAltAlleles
	 *            the targeted max alternative alleles.
	 * @return {@code true} iff this calculator implementation satisfies both
	 *         requirements.
	 */
	public boolean usableForParams(final int requestedPloidy, final int requestedMaxAltAlleles) {
		return (requiredPloidy == UNBOUND_PLOIDY || requiredPloidy == requestedPloidy)
				&& (maxAltAlleles == UNBOUND_ALTERNATIVE_ALLELE_COUNT || maxAltAlleles >= requestedMaxAltAlleles);
	}

	/**
	 * Resolve the constructor to use to instantiate calculators.
	 *
	 * @param clazz
	 *            target class. Assume not to be {@code null}.
	 */
	private Constructor<? extends AFCalculator> findInstantiationConstructor(
			final Class<? extends AFCalculator> clazz) {
		if (Modifier.isAbstract(clazz.getModifiers()))
			throw new IllegalStateException("AF calculator implementation class cannot be abstract");

		final Constructor<? extends AFCalculator> result;
		try {
			result = clazz.getDeclaredConstructor();
		} catch (final NoSuchMethodException e) {
			throw new IllegalStateException(
					"cannot find a suitable (int,int) constructor for the AFCalculator implementation " + this
							+ " class " + clazz.getName());
		}

		// Check whether there will be issue calling the constructor just due to
		// protections:
		if (Modifier.isPrivate(result.getModifiers())
				|| (!Modifier.isPublic(result.getModifiers()) && !clazz.getPackage().equals(getClass().getPackage())))
			throw new IllegalStateException("triple int constructor for AFCalculator implementation " + this + " class "
					+ clazz.getName() + " is not public ");
		return result;
	}

	/**
	 * Creates new instance.
	 *
	 * @throws IllegalStateException
	 *             if the instance could not be create due to some exception.
	 *             The {@link Exception#getCause() cause} will hold a reference
	 *             to the actual exception.
	 * @return never {@code null}.
	 */
	public AFCalculator  newInstance() {
		try {
			return constructor.newInstance();
		} catch (final Throwable e) {
			throw new IllegalStateException("could not instantiate AFCalculator for implementation " + this + " class "
					+ calculatorClass.getName());
		}
	}

	/**
	 * Returns the best (fastest) model give the required ploidy and alternative
	 * allele count.
	 *
	 * @param requiredPloidy
	 *            required ploidy
	 * @param requiredAlternativeAlleleCount
	 *            required alternative allele count.
	 * @param preferred
	 *            a preferred mode if any. A {@code null} indicate that we
	 *            should be try to use the default instead.
	 * @return never {@code null}
	 */
	public static AFCalculatorImplementation bestValue(final int requiredPloidy,
			final int requiredAlternativeAlleleCount, final AFCalculatorImplementation preferred) {
		final AFCalculatorImplementation preferredValue = preferred == null ? DEFAULT : preferred;
		if (preferredValue.usableForParams(requiredPloidy, requiredAlternativeAlleleCount))
			return preferredValue;
		else if (EXACT_INDEPENDENT.usableForParams(requiredPloidy, requiredAlternativeAlleleCount))
			return EXACT_INDEPENDENT;
		else if (EXACT_REFERENCE.usableForParams(requiredPloidy, requiredAlternativeAlleleCount))
			return EXACT_REFERENCE;
		else if (EXACT_GENERAL_INDEPENDENT.usableForParams(requiredPloidy, requiredAlternativeAlleleCount))
			return EXACT_GENERAL_INDEPENDENT;
		else
			return EXACT_GENERAL_PLOIDY;
	}

	/**
	 * Returns the value that corresponds to a given implementation calculator
	 * class.
	 *
	 * @param clazz
	 *            the target class.
	 *
	 * @throws IllegalArgumentException
	 *             if {@code clazz} is {@code null} or if it is abstract.
	 * @throws IllegalStateException
	 *             if
	 *
	 * @return never {@code null}.
	 */
	public static AFCalculatorImplementation fromCalculatorClass(final Class<? extends AFCalculator> clazz) {
		if (clazz == null)
			throw new IllegalArgumentException("input class cannot be null");
		final AFCalculatorImplementation result = calculatorClassToValue.get(clazz);
		if (result == null)
			throw new IllegalStateException(
					"Attempt to retrieve AFCalculatorImplementation instance from a non-registered calculator class "
							+ clazz.getName());
		return result;
	}

	// Initializes the content of the class to value map.
	private static Map<Class<? extends AFCalculator>, AFCalculatorImplementation> buildCalculatorClassToValueMap() {
		final Map<Class<? extends AFCalculator>, AFCalculatorImplementation> result = new HashMap<>(values().length);
		for (final AFCalculatorImplementation value : values())
			if (result.put(value.calculatorClass, value) != null)
				throw new IllegalStateException(
						"more than one value associated with class " + value.calculatorClass.getName());
		return result;
	}
}
