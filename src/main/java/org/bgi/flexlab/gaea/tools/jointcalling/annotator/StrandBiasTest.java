package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;

public abstract class StrandBiasTest extends InfoFieldAnnotation implements ActiveRegionBasedAnnotation {
	private static boolean getTableFromSamplesWarningLogged = false;
	private static boolean decodeSBBSWarningLogged = false;

	protected static final int ARRAY_DIM = 2;
	protected static final int ARRAY_SIZE = ARRAY_DIM * ARRAY_DIM;

	@Override
	public void initialize(Set<VCFHeaderLine> headerLines,Set<String> sampleList) {

	}

	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker, ChromosomeInformationShare ref, VariantContext vc) {

		// do not process if not a variant
		if (!vc.isVariant())
			return null;

		// if the genotype and strand bias are provided, calculate the
		// annotation from the Genotype (GT) field
		if (vc.hasGenotypes()) {
			for (final Genotype g : vc.getGenotypes()) {
				if (g.hasAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY)) {
					return calculateAnnotationFromGTfield(vc.getGenotypes());
				}
			}
		}

		return null;
	}

	protected abstract Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes);

	/**
	 */
	protected int[][] getTableFromSamples(final GenotypesContext genotypes, final int minCount) {
		if (genotypes == null) {
			if (!getTableFromSamplesWarningLogged) {
				getTableFromSamplesWarningLogged = true;
			}
			return null;
		}

		final int[] sbArray = { 0, 0, 0, 0 }; // reference-forward-reverse -by-
												// alternate-forward-reverse
		boolean foundData = false;

		for (final Genotype g : genotypes) {
			if (g.isNoCall() || !g.hasAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY))
				continue;

			foundData = true;
			int[] data;
			if (g.getAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY).getClass().equals(String.class)) {
				final String sbbsString = (String) g.getAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
				data = encodeSBBS(sbbsString);
			} else if (g.getAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY).getClass()
					.equals(ArrayList.class)) {
				ArrayList sbbsList = (ArrayList) g.getAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
				data = encodeSBBS(sbbsList);
			} else
				throw new IllegalArgumentException(
						"Unexpected " + GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY + " type");

			if (passesMinimumThreshold(data, minCount)) {
				for (int index = 0; index < sbArray.length; index++) {
					sbArray[index] += data[index];
				}
			}
		}

		return (foundData ? decodeSBBS(sbArray) : null);
	}

	/**
	 * Allocate and fill a 2x2 strand contingency table. In the end, it'll look
	 * something like this: fw rc allele1 # # allele2 # #
	 * 
	 * @return a 2x2 contingency table
	 */
	protected static int[][] getSNPContingencyTable(
			final Allele ref, final List<Allele> allAlts, final int minQScoreToConsider, final int minCount) {
		int[][] table = new int[ARRAY_DIM][ARRAY_DIM];

		return table;
	}

	/**
	 * Allocate and fill a 2x2 strand contingency table. In the end, it'll look
	 * something like this: fw rc allele1 # # allele2 # #
	 * 
	 * @return a 2x2 contingency table
	 */
	public static int[][] getContingencyTable(final VariantContext vc, final int minCount) {
		return null;
	}

	/**
	 * Helper method to copy the per-sample table to the main table
	 *
	 * @param perSampleTable
	 *            per-sample table (single dimension)
	 * @param mainTable
	 *            main table (two dimensions)
	 */
	private static void copyToMainTable(final int[] perSampleTable, final int[][] mainTable) {
		mainTable[0][0] += perSampleTable[0];
		mainTable[0][1] += perSampleTable[1];
		mainTable[1][0] += perSampleTable[2];
		mainTable[1][1] += perSampleTable[3];
	}

	/**
	 * Does this strand data array pass the minimum threshold for inclusion?
	 *
	 * @param data
	 *            the array
	 * @param minCount
	 *            The minimum threshold of counts in the array
	 * @return true if it passes the minimum threshold, false otherwise
	 */
	protected static boolean passesMinimumThreshold(final int[] data, final int minCount) {
		// the ref and alt totals must be greater than MIN_COUNT
		return data[0] + data[1] + data[2] + data[3] > minCount;
	}

	/**
	 * Helper function to parse the genotype annotation into the SB annotation
	 * array
	 * 
	 * @param string
	 *            the string that is returned by genotype.getAnnotation("SB")
	 * @return the array used by the per-sample Strand Bias annotation
	 */
	private static int[] encodeSBBS(final String string) {
		final int[] array = new int[ARRAY_SIZE];
		final StringTokenizer tokenizer = new StringTokenizer(string, ",", false);
		for (int index = 0; index < ARRAY_SIZE; index++) {
			array[index] = Integer.parseInt(tokenizer.nextToken());
		}
		return array;
	}

	/**
	 * Helper function to parse the genotype annotation into the SB annotation
	 * array
	 * 
	 * @param arrayList
	 *            the ArrayList returned from StrandBiasBySample.annotate()
	 * @return the array used by the per-sample Strand Bias annotation
	 */
	private static int[] encodeSBBS(final ArrayList<Integer> arrayList) {
		final int[] array = new int[ARRAY_SIZE];
		int index = 0;
		for (Integer item : arrayList)
			array[index++] = item.intValue();

		return array;
	}

	/**
	 * Helper function to turn the SB annotation array into a contingency table
	 * 
	 * @param array
	 *            the array used by the per-sample Strand Bias annotation
	 * @return the table used by the StrandOddsRatio annotation
	 */
	private static int[][] decodeSBBS(final int[] array) {
		if (array.length != ARRAY_SIZE) {
			if (!decodeSBBSWarningLogged) {
				decodeSBBSWarningLogged = true;
			}
			return null;
		}
		final int[][] table = new int[ARRAY_DIM][ARRAY_DIM];
		table[0][0] = array[0];
		table[0][1] = array[1];
		table[1][0] = array[2];
		table[1][1] = array[3];
		return table;
	}
}
