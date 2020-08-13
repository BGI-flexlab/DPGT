package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.genotyer.genotypecaller.ExactACcounts;
import org.bgi.flexlab.gaea.tools.genotyer.genotypecaller.ExactACset;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeAlleleCounts;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculator;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculators;
import org.bgi.flexlab.gaea.util.MathUtils;

import htsjdk.variant.variantcontext.Allele;

public class GeneralPloidyGenotypeLikelihoods {
	private static final int MAX_NUM_ALLELES_TO_CACHE = 20;
	private static final int MAX_NUM_SAMPLES_PER_POOL = 1000;

	private final static int[][] GenotypeLikelihoodVectorSizes = fillGLVectorSizeCache(MAX_NUM_ALLELES_TO_CACHE,
			2 * MAX_NUM_SAMPLES_PER_POOL);

	private static int[][] fillGLVectorSizeCache(int maxAlleles, int maxPloidy) {

		int[][] cache = new int[maxAlleles][maxPloidy];
		for (int numAlleles = 1; numAlleles < maxAlleles; numAlleles++) {
			for (int ploidy = 0; ploidy < maxPloidy; ploidy++) {

				if (numAlleles == 1)
					cache[numAlleles][ploidy] = 1;
				else if (ploidy == 1)
					cache[numAlleles][ploidy] = numAlleles;
				else {
					int acc = 0;
					for (int k = 0; k <= ploidy; k++)
						acc += cache[numAlleles - 1][ploidy - k];

					cache[numAlleles][ploidy] = acc;
				}
			}
		}
		return cache;
	}

	// Static methods
	public static void updateACset(final int[] newSetCounts, final LinkedList<ExactACset> ACqueue,
			final HashMap<ExactACcounts, ExactACset> indexesToACset) {

		final ExactACcounts index = new ExactACcounts(newSetCounts);
		if (!indexesToACset.containsKey(index)) {
			ExactACset newSet = new ExactACset(1, index);
			indexesToACset.put(index, newSet);
			ACqueue.add(newSet);
		}
	}

	public static int[] getAlleleCountFromPLIndex(final int nAlleles, final int numChromosomes, final int PLindex) {

		final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(numChromosomes,
				nAlleles);
		final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(PLindex);
		return alleleCounts.alleleCountsByIndex(nAlleles - 1);
	}

	public static int getLinearIndex(int[] vectorIdx, int numAlleles, int ploidy) {

		if (ploidy <= 0)
			return 0;

		int linearIdx = 0;
		int cumSum = ploidy;
		for (int k = numAlleles - 1; k >= 1; k--) {
			int idx = vectorIdx[k];
			// how many blocks are before current position
			if (idx == 0)
				continue;
			for (int p = 0; p < idx; p++)
				linearIdx += getNumLikelihoodElements(k, cumSum - p);

			cumSum -= idx;
		}

		return linearIdx;

	}

	public static int getNumLikelihoodElements(int numAlleles, int ploidy) {
		return GenotypeLikelihoodVectorSizes[numAlleles][ploidy];
	}

	public static double[] subsetToAlleles(final double[] oldLikelihoods, final int numChromosomes,
			final List<Allele> originalAlleles, final List<Allele> allelesToSubset) {

		int newPLSize = GeneralPloidyGenotypeLikelihoods.getNumLikelihoodElements(allelesToSubset.size(),
				numChromosomes);
		double[] newPLs = new double[newPLSize];

		int idx = 0;
		// First fill boolean array stating whether each original allele is
		// present in new mapping
		final boolean[] allelePresent = new boolean[originalAlleles.size()];
		for (Allele allele : originalAlleles)
			allelePresent[idx++] = allelesToSubset.contains(allele);

		// compute mapping from old idx to new idx
		// This might be needed in case new allele set is not ordered in the
		// same way as old set
		// Example. Original alleles: {T*,C,G,A}. New alleles: {G,C}.
		// Permutation key = [2,1]

		int[] permutationKey = new int[allelesToSubset.size()];
		for (int k = 0; k < allelesToSubset.size(); k++)
			// for each allele to subset, find corresponding index in original
			// allele list
			permutationKey[k] = originalAlleles.indexOf(allelesToSubset.get(k));

		final SumIterator iterator = new SumIterator(originalAlleles.size(), numChromosomes);

		while (iterator.hasNext()) {
			// for each entry in logPL table, associated originally with allele
			// count stored in vec[],
			// see if this allele count conformation will be present in new
			// logPL table.
			// For entry to be present, elements in dimensions not present in
			// requested allele list have to have count = 0
			int[] pVec = iterator.getCurrentVector();
			double pl = oldLikelihoods[iterator.getLinearIndex()];

			boolean keyPresent = true;
			for (int k = 0; k < allelePresent.length; k++)
				if (pVec[k] > 0 && !allelePresent[k])
					keyPresent = false;

			if (keyPresent) {// skip to next entry in logPLs if this
								// conformation is not present in subset

				final int[] newCount = new int[allelesToSubset.size()];

				// map from old allele mapping count to new allele mapping
				// In pseudo-Matlab notation: newCount = vec[permutationKey] for
				// permutationKey vector
				for (idx = 0; idx < newCount.length; idx++)
					newCount[idx] = pVec[permutationKey[idx]];

				// get corresponding index from new count
				int outputIdx = GeneralPloidyGenotypeLikelihoods.getLinearIndex(newCount, allelesToSubset.size(),
						numChromosomes);
				newPLs[outputIdx] = pl;
			}
			iterator.next();
		}

		return newPLs;
	}

	/**
	 * Crucial inner class that handles addressing elements of pool likelihoods.
	 * We store likelihoods as a map of form int[] -> double (to be more
	 * precise, IntArrayWrapper -> Double). For a given ploidy (chromosome
	 * count) and number of alleles, we need a form to iterate deterministically
	 * across all possible allele conformations. Problem equivalent to listing
	 * in deterministic order all possible ways in which N integers will sum to
	 * P, where N is number of alleles and P is number of chromosomes. There's
	 * an option to list all integers so that sum will be UP to P. For example,
	 * with P=2,N=2, restrictSumTo = 2 iterator will produce [2 0] [1 1] [0 2]
	 *
	 *
	 */
	public static class SumIterator {
		private int[] currentState;
		private final int[] finalState;
		private final int restrictSumTo;
		private final int dim;
		private boolean hasNext;
		private int linearIndex;
		private int currentSum;

		/**
		 * Default constructor. Typical use case: restrictSumTo = -1 if there's
		 * no sum restriction, or will generate int[] vectors so that all add to
		 * this value.
		 *
		 * @param finalState
		 *            End state - typically we should set value to (P,P,P,...)
		 * @param restrictSumTo
		 *            See above
		 */
		public SumIterator(final int[] finalState, final int restrictSumTo) {
			this.finalState = finalState;
			this.dim = finalState.length;
			this.restrictSumTo = restrictSumTo;
			currentState = new int[dim];
			reset();

		}

		/**
		 * Shortcut constructor for common use case: iterator will produce all
		 * vectors of length numAlleles whose sum = numChromosomes
		 * 
		 * @param numAlleles
		 *            Number of alleles
		 * @param numChromosomes
		 *            Ploidy
		 */
		public SumIterator(final int numAlleles, final int numChromosomes) {
			this(getInitialStateVector(numAlleles, numChromosomes), numChromosomes);
		}

		private static int[] getInitialStateVector(final int nAlleles, final int numChromosomes) {
			int[] initialState = new int[nAlleles];
			Arrays.fill(initialState, numChromosomes);
			return initialState;
		}

		public void setInitialStateVector(final int[] stateVector) {
			if (restrictSumTo > 0) {
				// check that desired vector is valid
				if (MathUtils.sum(stateVector) != restrictSumTo)
					throw new UserException("BUG: initial state vector nor compatible with sum iterator");

				final int numAlleles = currentState.length;
				final int ploidy = restrictSumTo;

				linearIndex = GeneralPloidyGenotypeLikelihoods.getLinearIndex(stateVector, numAlleles, ploidy);
			} else
				throw new UserException("BUG: Not supported");

		}

		public void next() {
			int initialDim = (restrictSumTo > 0) ? 1 : 0;
			hasNext = next(finalState, initialDim);
			if (hasNext)
				linearIndex++;
		}

		private boolean next(final int[] finalState, int initialDim) {
			boolean hasNextState = false;
			for (int currentDim = initialDim; currentDim < finalState.length; currentDim++) {
				final int x = currentState[currentDim] + 1;

				if (x > finalState[currentDim] || (currentSum >= restrictSumTo && initialDim > 0)) {
					// update vector sum, and reset position
					currentSum -= currentState[currentDim];
					currentState[currentDim] = 0;
					if (currentDim >= dim - 1) {
						hasNextState = false;
						break;
					}
				} else {
					currentState[currentDim] = x;
					hasNextState = true;
					currentSum++;
					break;
				}
			}
			if (initialDim > 0) {
				currentState[0] = restrictSumTo - currentSum;
			}
			return hasNextState;
		}

		public void reset() {
			Arrays.fill(currentState, 0);
			if (restrictSumTo > 0)
				currentState[0] = restrictSumTo;
			hasNext = true;
			linearIndex = 0;
			currentSum = 0;
		}

		public int[] getCurrentVector() {
			return currentState;
		}

		public int[] getCurrentAltVector() {
			return Arrays.copyOfRange(currentState, 1, currentState.length);
		}

		/*
		 * public int getCurrentSum() { return currentSum; }
		 */
		public int getLinearIndex() {
			return linearIndex;
		}

		public boolean hasNext() {
			return hasNext;
		}
	}
}
