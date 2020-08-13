package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

import org.bgi.flexlab.gaea.tools.jointcalling.util.SimpleTimer;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public abstract class AFCalculator implements Cloneable {
	public static final int MAX_NUM_PL_VALUES_DEFAULT = 100;

	protected int maxNumPLValues = MAX_NUM_PL_VALUES_DEFAULT; // if PL vectors longer than this # of elements,
																// don't log them
	protected static int maxNumPLValuesObserved = 0;
	protected static long numTimesMaxNumPLValuesExceeded = 0;

	private SimpleTimer callTimer = new SimpleTimer();
	private StateTracker stateTracker;
	private ExactCallLogger exactCallLogger = null;

	protected AFCalculator() {
	}

	/**
	 * Enable exact call logging to file
	 *
	 * @param exactCallsLog
	 *            the destination file
	 */
	public void enableProcessLog(final File exactCallsLog) {
		exactCallLogger = new ExactCallLogger(exactCallsLog);
	}

	/**
	 * Set the maximum number of PL values to log. If the number of PL values
	 * exceeds this, no PL values will be logged.
	 * 
	 * @param maxNumPLValues
	 *            maximum number of PL values to log
	 */
	public AFCalculator setMaxNumPLValues(final int maxNumPLValues) {
		this.maxNumPLValues = maxNumPLValues;
		return this;
	}

	/**
	 * Compute the probability of the alleles segregating given the genotype
	 * likelihoods of the samples in vc
	 *
	 * @param vc
	 *            the VariantContext holding the alleles and sample information.
	 *            The VariantContext must have at least 1 alternative allele
	 * @param log10AlleleFrequencyPriors
	 *            a prior vector nSamples x 2 in length indicating the Pr(AF =
	 *            i)
	 * @return result (for programming convenience)
	 */
	public AFCalculationResult getLog10PNonRef(final VariantContext vc, final int defaultPloidy,
			final int maximumAlternativeAlleles, final double[] log10AlleleFrequencyPriors) {
		SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
		if (vc == null)
			throw new IllegalArgumentException("VariantContext cannot be null");
		if (vc.getNAlleles() == 1)
			throw new IllegalArgumentException(
					"VariantContext has only a single reference allele, but getLog10PNonRef requires at least one at all "
							+ vc);
		if (log10AlleleFrequencyPriors == null)
			throw new IllegalArgumentException("priors vector cannot be null");

		// reset the result, so we can store our new result there
		final StateTracker stateTracker = getStateTracker(true, maximumAlternativeAlleles);

		// TODO All implementations of the reduce-scope seems to employ a bad
		// criterion to
		// TODO decide what alleles to keep. This must be changed eventually.
		// TODO issue {@see
		// https://github.com/broadinstitute/gsa-unstable/issues/1376}
//		System.out.println(formatter.format(new Date())+"\tbefore reduceScope");
		final VariantContext vcWorking = reduceScope(vc, defaultPloidy, maximumAlternativeAlleles);
//		System.out.println(formatter.format(new Date())+"\tafter reduceScope, before computeLog10PNonRef");
		callTimer.start();
		final AFCalculationResult result = computeLog10PNonRef(vcWorking, defaultPloidy, log10AlleleFrequencyPriors,
				stateTracker);
//		System.out.println(formatter.format(new Date())+"\tafter computeLog10PNonRef");
		final long nanoTime = callTimer.getElapsedTimeNano();

		if (exactCallLogger != null)
			exactCallLogger.printCallInfo(vcWorking, log10AlleleFrequencyPriors, nanoTime, result);
//		System.out.println(formatter.format(new Date())+"\tafter printCallInfo");
		return result;
	}

	/**
	 * Convert the final state of the state tracker into our result as an
	 * AFCalculationResult
	 *
	 * Assumes that stateTracker has been updated accordingly
	 *
	 * @param vcWorking
	 *            the VariantContext we actually used as input to the calc model
	 *            (after reduction)
	 * @param log10AlleleFrequencyPriors
	 *            the priors by AC vector
	 * @return a AFCalculationResult describing the result of this calculation
	 */
	protected AFCalculationResult getResultFromFinalState(final VariantContext vcWorking,
			final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
		stateTracker.setAllelesUsedInGenotyping(vcWorking.getAlleles());
		return stateTracker.toAFCalculationResult(log10AlleleFrequencyPriors);
	}

	// ---------------------------------------------------------------------------
	//
	// Abstract methods that should be implemented by concrete implementations
	// to actually calculate the AF
	//
	// ---------------------------------------------------------------------------

	/**
	 * Look at VC and perhaps return a new one of reduced complexity, if that's
	 * necessary
	 *
	 * Used before the call to computeLog10PNonRef to simply the calculation job
	 * at hand, if vc exceeds bounds. For example, if VC has 100 alt alleles
	 * this function may decide to only genotype the best 2 of them.
	 *
	 * @param vc
	 *            the initial VC provided by the caller to this AFcalculation
	 * @return a potentially simpler VC that's more tractable to genotype
	 */
	protected abstract VariantContext reduceScope(final VariantContext vc, final int defaultPloidy,
			final int maximumAlternativeAlleles);

	/**
	 * Actually carry out the log10PNonRef calculation on vc, storing results in
	 * results
	 *
	 * @param vc
	 *            variant context with alleles and genotype likelihoods, must
	 *            have at least one alt allele
	 * @param log10AlleleFrequencyPriors
	 *            priors
	 * @return a AFCalcResult object describing the results of this calculation
	 */
	protected abstract AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
			final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker);

	/**
	 * Subset VC to the just allelesToUse, updating genotype likelihoods
	 *
	 * Must be overridden by concrete subclasses
	 *
	 * @param vc
	 *            variant context with alleles and genotype likelihoods
	 * @param defaultPloidy
	 *            default ploidy to assume in case {@code vc} does not indicate
	 *            it for a sample.
	 * @param allelesToUse
	 *            alleles to subset
	 * @param assignGenotypes
	 * @return GenotypesContext object
	 */
	public abstract GenotypesContext subsetAlleles(final VariantContext vc, final int defaultPloidy,
			final List<Allele> allelesToUse, final boolean assignGenotypes);

	// ---------------------------------------------------------------------------
	//
	// accessors
	//
	// ---------------------------------------------------------------------------

	/**
	 * Retrieves the state tracker.
	 *
	 * <p>
	 * The tracker will be reset if so requested or if it needs to be resized
	 * due to an increase in the maximum number of alleles is must be able to
	 * handle.
	 * </p>
	 *
	 * @param reset
	 *            make sure the tracker is reset.
	 * @param maximumAlternativeAlleleCount
	 *            the maximum alternative allele count it must be able to
	 *            handle. Has no effect if the current tracker is able to handle
	 *            that number.
	 *
	 * @return never {@code null}
	 */
	protected StateTracker getStateTracker(final boolean reset, final int maximumAlternativeAlleleCount) {
		if (stateTracker == null)
			stateTracker = new StateTracker(maximumAlternativeAlleleCount);
		else if (reset)
			stateTracker.reset(maximumAlternativeAlleleCount);
		else
			stateTracker.ensureMaximumAlleleCapacity(maximumAlternativeAlleleCount);
		return stateTracker;
	}

	/**
	 * Used by testing code.
	 *
	 * Please don't use this method in production.
	 *
	 * @deprecated
	 */
	@Deprecated
	protected int getAltAlleleCountOfMAP(final int allele) {
		return getStateTracker(false, allele + 1).getAlleleCountsOfMAP()[allele];
	}

	/**
	 * Strips PLs from the specified GenotypeBuilder if their number exceeds the
	 * maximum allowed. Corresponding counters are updated.
	 * 
	 * @param gb
	 *            the GenotypeBuilder to modify
	 * @param vc
	 *            the VariantContext
	 * @param sampleName
	 *            the sample name
	 * @param newLikelihoods
	 *            the PL array
	 */
	protected void removePLsIfMaxNumPLValuesExceeded(final GenotypeBuilder gb, final VariantContext vc,
			final String sampleName, final double[] newLikelihoods) {
		final int numPLValuesFound = newLikelihoods.length;
		if (numPLValuesFound > maxNumPLValues) {
			numTimesMaxNumPLValuesExceeded++;
			gb.noPL();
			if (numPLValuesFound > maxNumPLValuesObserved) {
				maxNumPLValuesObserved = numPLValuesFound;
			}
		}
	}

	/**
	 * Logs the number of times the maximum allowed number of PLs was exceeded
	 * and the largest number of PLs observed. The corresponding counters are
	 * reset.
	 */
	public void printFinalMaxNumPLValuesWarning() {
		maxNumPLValuesObserved = 0;
		numTimesMaxNumPLValuesExceeded = 0;
	}
}
