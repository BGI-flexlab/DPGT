package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PairHMMIndelErrorModel;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;

public class AnnotationUtils {
	public static final String ANNOTATION_HC_WARN_MSG = " annotation will not be calculated, must be called from HaplotypeCaller or MuTect2";
	public static final int WARNINGS_LOGGED_SIZE = 3;

	/**
	 * Helper function to convert a List of Doubles to a comma-separated String
	 * 
	 * @param valueList
	 *            the ArrayList with Double data
	 * @return a comma-separated String
	 */
	public static String encodeValueList(final List<Double> valueList, final String precisionFormat) {
		List<String> outputList = new ArrayList<>();
		for (Double d : valueList) {
			outputList.add(String.format(precisionFormat, d));
		}
		return StringUtils.join(outputList, ",");
	}

	/**
	 * Helper function to convert a List of Strings to a comma-separated String
	 * 
	 * @param stringList
	 *            the ArrayList with String data
	 * @return a comma-separated String
	 */
	public static String encodeStringList(final List<String> stringList) {
		return StringUtils.join(stringList, ",");
	}

	// this method is intended to reconcile uniquified sample names
	// it comes into play when calling this annotation from GenotypeGVCFs with
	// --uniquifySamples because founderIds
	// is derived from the sampleDB, which comes from the input sample names,
	// but vc will have uniquified (i.e. different)
	// sample names. Without this check, the founderIds won't be found in the vc
	// and the annotation won't be calculated.
	public static Set<String> validateFounderIDs(final Set<String> founderIds, final VariantContext vc) {
		Set<String> vcSamples = new HashSet<>();
		Set<String> returnIDs = founderIds;
		vcSamples.addAll(vc.getSampleNames());
		if (!vcSamples.isEmpty()) {
			if (founderIds != null) {
				vcSamples.removeAll(founderIds);
				if (vcSamples.equals(vc.getSampleNames()))
					returnIDs = vc.getSampleNames();
			}
		}
		return returnIDs;
	}

	/**
	 * Get the position of a variant within a read with respect to the closer
	 * end, accounting for hard clipped bases and low quality ends Used by
	 * ReadPosRankSum annotations
	 *
	 * @param read
	 *            a read containing the variant
	 * @param initialReadPosition
	 *            the position based on the modified, post-hard-clipped CIGAR
	 * @return read position
	 */
	public static int getFinalVariantReadPosition(final SAMRecord read, final int initialReadPosition) {
		final int numAlignedBases = getNumAlignedBases(read);

		int readPos = initialReadPosition;
		// TODO: this doesn't work for the middle-right position if we index
		// from zero
		if (initialReadPosition > numAlignedBases / 2) {
			readPos = numAlignedBases - (initialReadPosition + 1);
		}
		return readPos;

	}

	/**
	 *
	 * @param read
	 *            a read containing the variant
	 * @return the number of hard clipped and low qual bases at the read start
	 *         (where start is the leftmost end w.r.t. the reference)
	 */
	public static int getNumClippedBasesAtStart(final SAMRecord read) {
		// check for hard clips (never consider these bases):
		final Cigar c = read.getCigar();
		final CigarElement first = c.getCigarElement(0);

		int numStartClippedBases = 0;
		if (first.getOperator() == CigarOperator.H) {
			numStartClippedBases = first.getLength();
		}
		final byte[] unclippedReadBases = read.getReadBases();
		final byte[] unclippedReadQuals = read.getBaseQualities();

		// Do a stricter base clipping than provided by CIGAR string, since this
		// one may be too conservative,
		// and may leave a string of Q2 bases still hanging off the reads.
		// TODO: this code may not even get used because HaplotypeCaller already
		// hard clips low quality tails
		for (int i = numStartClippedBases; i < unclippedReadBases.length; i++) {
			if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
				numStartClippedBases++;
			else
				break;

		}

		return numStartClippedBases;
	}

	/**
	 *
	 * @param read
	 *            a read containing the variant
	 * @return number of non-hard clipped, aligned bases (excluding low quality
	 *         bases at either end)
	 */
	// TODO: this is bizarre -- this code counts hard clips, but then subtracts
	// them from the read length, which already doesn't count hard clips
	public static int getNumAlignedBases(final SAMRecord read) {
		return read.getReadLength() - getNumClippedBasesAtStart(read) - getNumClippedBasesAtEnd(read);
	}

	/**
	 *
	 * @param read
	 *            a read containing the variant
	 * @return number of hard clipped and low qual bases at the read end (where
	 *         end is right end w.r.t. the reference)
	 */
	public static int getNumClippedBasesAtEnd(final SAMRecord read) {
		// check for hard clips (never consider these bases):
		final Cigar c = read.getCigar();
		CigarElement last = c.getCigarElement(c.numCigarElements() - 1);

		int numEndClippedBases = 0;
		if (last.getOperator() == CigarOperator.H) {
			numEndClippedBases = last.getLength();
		}
		final byte[] unclippedReadBases = read.getReadBases();
		final byte[] unclippedReadQuals = read.getBaseQualities();

		// Do a stricter base clipping than provided by CIGAR string, since this
		// one may be too conservative,
		// and may leave a string of Q2 bases still hanging off the reads.
		// TODO: this code may not even get used because HaplotypeCaller already
		// hard clips low quality tails
		for (int i = unclippedReadBases.length - numEndClippedBases - 1; i >= 0; i--) {
			if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
				numEndClippedBases++;
			else
				break;
		}

		return numEndClippedBases;
	}
}
