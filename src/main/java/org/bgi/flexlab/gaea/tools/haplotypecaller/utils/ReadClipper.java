package org.bgi.flexlab.gaea.tools.haplotypecaller.utils;

import java.util.ArrayList;
import java.util.List;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.ReadUtils;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class ReadClipper {
	final GaeaSamRecord read;
	boolean wasClipped;
	List<ClippingOp> ops = null;

	/**
	 * Initializes a ReadClipper object.
	 *
	 * You can set up your clipping operations using the addOp method. When
	 * you're ready to generate a new read with all the clipping operations, use
	 * clipRead().
	 *
	 * Note: Use this if you want to set up multiple operations on the read
	 * using the ClippingOp class. If you just want to apply one of the typical
	 * modes of clipping, use the static clipping functions available in this
	 * class instead.
	 *
	 * @param read
	 *            the read to clip
	 */
	public ReadClipper(final GaeaSamRecord read) {
		Utils.nonNull(read);
		this.read = read;
		this.wasClipped = false;
	}

	/**
	 * Add clipping operation to the read.
	 *
	 * You can add as many operations as necessary to this read before clipping.
	 * Beware that the order in which you add these operations matter. For
	 * example, if you hard clip the beginning of a read first then try to hard
	 * clip the end, the indices will have changed. Make sure you know what
	 * you're doing, otherwise just use the static functions below that take
	 * care of the ordering for you.
	 *
	 * Note: You only choose the clipping mode when you use clipRead()
	 *
	 * @param op
	 *            a ClippingOp object describing the area you want to clip.
	 */
	public void addOp(final ClippingOp op) {
		Utils.nonNull(op);
		if (ops == null) {
			ops = new ArrayList<>();
		}
		ops.add(op);
	}

	/**
	 * Check the list of operations set up for this read.
	 *
	 * @return a list of the operations set up for this read.
	 */
	public List<ClippingOp> getOps() {
		return ops;
	}

	/**
	 * Check whether or not this read has been clipped.
	 * 
	 * @return true if this read has produced a clipped read, false otherwise.
	 */
	public boolean wasClipped() {
		return wasClipped;
	}

	/**
	 * The original read.
	 *
	 * @return returns the read to be clipped (original)
	 */
	public GaeaSamRecord getRead() {
		return read;
	}

	/**
	 * Clips a read according to ops and the chosen algorithm.
	 *
	 * @param algorithm
	 *            What mode of clipping do you want to apply for the stacked
	 *            operations.
	 * @return the read with the clipping applied.
	 */
	public GaeaSamRecord clipRead(final ClippingRepresentation algorithm) {
		return clipRead(algorithm, true);
	}

	private GaeaSamRecord clipRead(final ClippingRepresentation algorithm, boolean runAsserts) {
		Utils.nonNull(algorithm);
		if (ops == null) {
			return getRead();
		}

		GaeaSamRecord clippedRead = read;
		for (final ClippingOp op : getOps()) {
			final int readLength = clippedRead.getReadLength();
			// check if the clipped read can still be clipped in the range
			// requested
			if (op.start < readLength) {
				ClippingOp fixedOperation = op;
				if (op.stop >= readLength) {
					fixedOperation = new ClippingOp(op.start, readLength - 1);
				}

				clippedRead = fixedOperation.apply(algorithm, clippedRead, runAsserts);
			}
		}
		wasClipped = true;
		ops.clear();
		if (clippedRead.isEmpty()) {
			return ReadUtils.emptyRead(clippedRead);
		}
		return clippedRead;
	}

	/**
	 * Hard clips the left tail of a read up to (and including) refStop using
	 * reference coordinates.
	 *
	 * @param refStop
	 *            the last base to be hard clipped in the left tail of the read.
	 * @return a new read, without the left tail.
	 */
	private GaeaSamRecord hardClipByReferenceCoordinatesLeftTail(final int refStop) {
		return clipByReferenceCoordinates(-1, refStop, ClippingRepresentation.HARDCLIP_BASES, true);
	}

	public static GaeaSamRecord hardClipByReferenceCoordinatesLeftTail(final GaeaSamRecord read, final int refStop) {
		return (new ReadClipper(read)).clipByReferenceCoordinates(-1, refStop, ClippingRepresentation.HARDCLIP_BASES,
				true);
	}

	/**
	 * Hard clips the right tail of a read starting at (and including) refStart
	 * using reference coordinates.
	 *
	 * @param refStart
	 *            refStop the first base to be hard clipped in the right tail of
	 *            the read.
	 * @return a new read, without the right tail.
	 */
	private GaeaSamRecord hardClipByReferenceCoordinatesRightTail(final int refStart) {
		return clipByReferenceCoordinates(refStart, -1, ClippingRepresentation.HARDCLIP_BASES, true);
	}

	public static GaeaSamRecord hardClipByReferenceCoordinatesRightTail(final GaeaSamRecord read, final int refStart) {
		return (new ReadClipper(read)).clipByReferenceCoordinates(refStart, -1, ClippingRepresentation.HARDCLIP_BASES,
				true);
	}

	/**
	 * Hard clips a read using read coordinates.
	 *
	 * @param start
	 *            the first base to clip (inclusive)
	 * @param stop
	 *            the last base to clip (inclusive)
	 * @return a new read, without the clipped bases
	 */
	private GaeaSamRecord hardClipByReadCoordinates(final int start, final int stop) {
		if (read.isEmpty() || (start == 0 && stop == read.getReadLength() - 1)) {
			return ReadUtils.emptyRead(read);
		}

		this.addOp(new ClippingOp(start, stop));
		return clipRead(ClippingRepresentation.HARDCLIP_BASES);
	}

	public static GaeaSamRecord hardClipByReadCoordinates(final GaeaSamRecord read, final int start, final int stop) {
		return (new ReadClipper(read)).hardClipByReadCoordinates(start, stop);
	}

	/**
	 * Hard clips both tails of a read. Left tail goes from the beginning to the
	 * 'left' coordinate (inclusive) Right tail goes from the 'right' coordinate
	 * (inclusive) until the end of the read
	 *
	 * @param left
	 *            the coordinate of the last base to be clipped in the left tail
	 *            (inclusive)
	 * @param right
	 *            the coordinate of the first base to be clipped in the right
	 *            tail (inclusive)
	 * @return a new read, without the clipped bases
	 */
	private GaeaSamRecord hardClipBothEndsByReferenceCoordinates(final int left, final int right) {
		if (read.isEmpty() || left == right) {
			return ReadUtils.emptyRead(read);
		}
		final GaeaSamRecord leftTailRead = clipByReferenceCoordinates(right, -1, ClippingRepresentation.HARDCLIP_BASES,
				true);

		// after clipping one tail, it is possible that the consequent hard
		// clipping of adjacent deletions
		// make the left cut index no longer part of the read. In that case,
		// clip the read entirely.
		if (left > leftTailRead.getEnd()) {
			return ReadUtils.emptyRead(read);
		}

		final ReadClipper clipper = new ReadClipper(leftTailRead);
		return clipper.hardClipByReferenceCoordinatesLeftTail(left);
	}

	public static GaeaSamRecord hardClipBothEndsByReferenceCoordinates(final GaeaSamRecord read, final int left,
			final int right) {
		return (new ReadClipper(read)).hardClipBothEndsByReferenceCoordinates(left, right);
	}

	/**
	 * Clips any contiguous tail (left, right or both) with base quality lower
	 * than lowQual using the desired algorithm.
	 *
	 * This function will look for low quality tails and hard clip them away. A
	 * low quality tail ends when a base has base quality greater than lowQual.
	 *
	 * @param algorithm
	 *            the algorithm to use (HardClip, SoftClip, Write N's,...)
	 * @param lowQual
	 *            every base quality lower than or equal to this in the tail of
	 *            the read will be hard clipped
	 * @return a new read without low quality tails
	 */
	private GaeaSamRecord clipLowQualEnds(final ClippingRepresentation algorithm, final byte lowQual) {
		if (read.isEmpty()) {
			return read;
		}

		final int readLength = read.getReadLength();
		int leftClipIndex = 0;
		int rightClipIndex = readLength - 1;

		// check how far we can clip both sides
		while (rightClipIndex >= 0 && read.getBaseQuality(rightClipIndex) <= lowQual) {
			rightClipIndex--;
		}
		while (leftClipIndex < readLength && read.getBaseQuality(leftClipIndex) <= lowQual) {
			leftClipIndex++;
		}

		// if the entire read should be clipped, then return an empty read.
		if (leftClipIndex > rightClipIndex) {
			return ReadUtils.emptyRead(read);
		}

		if (rightClipIndex < readLength - 1) {
			this.addOp(new ClippingOp(rightClipIndex + 1, readLength - 1));
		}
		if (leftClipIndex > 0) {
			this.addOp(new ClippingOp(0, leftClipIndex - 1));
		}
		return this.clipRead(algorithm);
	}

	private GaeaSamRecord hardClipLowQualEnds(final byte lowQual) {
		return this.clipLowQualEnds(ClippingRepresentation.HARDCLIP_BASES, lowQual);
	}

	public static GaeaSamRecord clipLowQualEnds(final GaeaSamRecord read, final byte lowQual,
			final ClippingRepresentation algorithm) {
		return (new ReadClipper(read)).clipLowQualEnds(algorithm, lowQual);
	}

	public static GaeaSamRecord hardClipLowQualEnds(final GaeaSamRecord read, final byte lowQual) {
		return (new ReadClipper(read)).hardClipLowQualEnds(lowQual);
	}

	/**
	 * Will hard clip every soft clipped bases in the read.
	 *
	 * @return a new read without the soft clipped bases
	 */
	private GaeaSamRecord hardClipSoftClippedBases() {
		if (read.isEmpty()) {
			return read;
		}

		int readIndex = 0;
		int cutLeft = -1; // first position to hard clip (inclusive)
		int cutRight = -1; // first position to hard clip (inclusive)
		boolean rightTail = false; // trigger to stop clipping the left tail and
									// start cutting the right tail

		for (final CigarElement cigarElement : read.getCigar().getCigarElements()) {
			if (cigarElement.getOperator() == CigarOperator.SOFT_CLIP) {
				if (rightTail) {
					cutRight = readIndex;
				} else {
					cutLeft = readIndex + cigarElement.getLength() - 1;
				}
			} else if (cigarElement.getOperator() != CigarOperator.HARD_CLIP) {
				rightTail = true;
			}

			if (cigarElement.getOperator().consumesReadBases()) {
				readIndex += cigarElement.getLength();
			}
		}

		// It is extremely important that we cut the end first otherwise the
		// read coordinates change.
		if (cutRight >= 0) {
			this.addOp(new ClippingOp(cutRight, read.getReadLength() - 1));
		}
		if (cutLeft >= 0) {
			this.addOp(new ClippingOp(0, cutLeft));
		}

		return clipRead(ClippingRepresentation.HARDCLIP_BASES);
	}

	public static GaeaSamRecord hardClipSoftClippedBases(final GaeaSamRecord read) {
		return new ReadClipper(read).hardClipSoftClippedBases();
	}

	/**
	 * Hard clip the read to the variable region (from refStart to refStop)
	 *
	 * @param read
	 *            the read to be clipped
	 * @param refStart
	 *            the beginning of the variant region (inclusive)
	 * @param refStop
	 *            the end of the variant region (inclusive)
	 * @return the read hard clipped to the variant region
	 */
	public static GaeaSamRecord hardClipToRegion(final GaeaSamRecord read, final int refStart, final int refStop) {
		final int start = read.getStart();
		final int stop = read.getEnd();
		return hardClipToRegion(read, refStart, refStop, start, stop);
	}

	/**
	 * Hard clip the read to the variable region (from refStart to refStop)
	 * processing also the clipped bases
	 *
	 * @param read
	 *            the read to be clipped
	 * @param refStart
	 *            the beginning of the variant region (inclusive)
	 * @param refStop
	 *            the end of the variant region (inclusive)
	 * @return the read hard clipped to the variant region
	 */
	public static GaeaSamRecord hardClipToRegionIncludingClippedBases(final GaeaSamRecord read, final int refStart,
			final int refStop) {
		final int start = read.getUnclippedStart();
		final int stop = start
				+ CigarUtils.countRefBasesBasedOnCigar(read, 0, read.getCigar().getCigarElements().size()) - 1;
		return hardClipToRegion(read, refStart, refStop, start, stop);
	}

	private static GaeaSamRecord hardClipToRegion(final GaeaSamRecord read, final int refStart, final int refStop,
			final int alignmentStart, final int alignmentStop) {
		// check if the read is contained in region
		if (alignmentStart <= refStop && alignmentStop >= refStart) {
			if (alignmentStart < refStart && alignmentStop > refStop) {
				return hardClipBothEndsByReferenceCoordinates(read, refStart - 1, refStop + 1);
			} else if (alignmentStart < refStart) {
				return hardClipByReferenceCoordinatesLeftTail(read, refStart - 1);
			} else if (alignmentStop > refStop) {
				return hardClipByReferenceCoordinatesRightTail(read, refStop + 1);
			}
			return read;
		} else {
			return ReadUtils.emptyRead(read);
		}

	}

	/**
	 * Checks if a read contains adaptor sequences. If it does, hard clips them
	 * out.
	 *
	 * Note: To see how a read is checked for adaptor sequence see
	 * ReadUtils.getAdaptorBoundary()
	 *
	 * @return a new read without adaptor sequence
	 */
	private GaeaSamRecord hardClipAdaptorSequence() {
		final int adaptorBoundary = ReadUtils.getAdaptorBoundary(read);

		if (adaptorBoundary == Integer.MIN_VALUE || !ReadUtils.isInsideRead(read, adaptorBoundary)) {
			return read;
		}

		return read.getReadNegativeStrandFlag() ? hardClipByReferenceCoordinatesLeftTail(adaptorBoundary)
				: hardClipByReferenceCoordinatesRightTail(adaptorBoundary);
	}

	public static GaeaSamRecord hardClipAdaptorSequence(final GaeaSamRecord read) {
		return new ReadClipper(read).hardClipAdaptorSequence();
	}

	/**
	 * Hard clips any leading insertions in the read. Only looks at the
	 * beginning of the read, not the end.
	 *
	 * @return a new read without leading insertions
	 */
	private GaeaSamRecord hardClipLeadingInsertions() {
		if (read.isEmpty()) {
			return read;
		}

		for (final CigarElement cigarElement : read.getCigar().getCigarElements()) {
			if (cigarElement.getOperator() != CigarOperator.HARD_CLIP
					&& cigarElement.getOperator() != CigarOperator.SOFT_CLIP
					&& cigarElement.getOperator() != CigarOperator.INSERTION) {
				break;
			} else if (cigarElement.getOperator() == CigarOperator.INSERTION) {
				this.addOp(new ClippingOp(0, cigarElement.getLength() - 1));
			}
		}
		return clipRead(ClippingRepresentation.HARDCLIP_BASES);
	}

	public static GaeaSamRecord hardClipLeadingInsertions(final GaeaSamRecord read) {
		return (new ReadClipper(read)).hardClipLeadingInsertions();
	}

	/**
	 * Turns soft clipped bases into matches
	 * 
	 * @return a new read with every soft clip turned into a match
	 */
	private GaeaSamRecord revertSoftClippedBases() {
		if (read.isEmpty()) {
			return read;
		}

		this.addOp(new ClippingOp(0, 0));
		return this.clipRead(ClippingRepresentation.REVERT_SOFTCLIPPED_BASES);
	}

	/**
	 * Reverts ALL soft-clipped bases
	 *
	 * @param read
	 *            the read
	 * @return the read with all soft-clipped bases turned into matches
	 */
	public static GaeaSamRecord revertSoftClippedBases(final GaeaSamRecord read) {
		return new ReadClipper(read).revertSoftClippedBases();
	}

	protected GaeaSamRecord hardClipByReferenceCoordinates(final int refStart, final int refStop) {
		return clipByReferenceCoordinates(refStart, refStop, ClippingRepresentation.HARDCLIP_BASES, true);
	}

	/**
	 * Generic functionality to clip a read, used internally by
	 * hardClipByReferenceCoordinatesLeftTail and
	 * hardClipByReferenceCoordinatesRightTail. Should not be used directly.
	 *
	 * Note, it REQUIRES you to give the directionality of your hard clip (i.e.
	 * whether you're clipping the left of right tail) by specifying either
	 * refStart < 0 or refStop < 0.
	 *
	 * @param refStart
	 *            first base to clip (inclusive)
	 * @param refStop
	 *            last base to clip (inclusive)
	 * @param clippingOp
	 *            clipping operation to be performed
	 * @return a new read, without the clipped bases
	 */
	protected GaeaSamRecord clipByReferenceCoordinates(final int refStart, final int refStop,
			ClippingRepresentation clippingOp, boolean runAsserts) {
		if (read.isEmpty()) {
			return read;
		}
		if ((clippingOp == ClippingRepresentation.SOFTCLIP_BASES) && read.isUnmapped()) {
			throw new UserException(
					"Cannot softclip read " + read.getReadName() + " by reference coordinates because it is unmapped");
		}

		final int start;
		final int stop;

		// Determine the read coordinate to start and stop hard clipping
		if (refStart < 0) {
			if (refStop < 0) {
				throw new UserException(
						"Only one of refStart or refStop must be < 0, not both (" + refStart + ", " + refStop + ")");
			}
			start = 0;
			stop = ReadUtils.getReadCoordinateForReferenceCoordinate(read, refStop, ReadUtils.ClippingTail.LEFT_TAIL);
		} else {
			if (refStop >= 0) {
				throw new UserException("Either refStart or refStop must be < 0 (" + refStart + ", " + refStop + ")");
			}
			start = ReadUtils.getReadCoordinateForReferenceCoordinate(read, refStart,
					ReadUtils.ClippingTail.RIGHT_TAIL);
			stop = read.getReadLength() - 1;
		}

		if (start < 0 || stop > read.getReadLength() - 1) {
			throw new UserException("Trying to clip before the start or after the end of a read");
		}

		if (start > stop) {
			throw new UserException(String.format(
					"START (%d) > (%d) STOP -- this should never happen, please check read: %s (CIGAR: %s)", start,
					stop, read, read.getCigar().toString()));
		}

		if (start > 0 && stop < read.getReadLength() - 1) {
			throw new UserException(String.format("Trying to clip the middle of the read: start %d, stop %d, cigar: %s",
					start, stop, read.getCigar().toString()));
		}
		this.addOp(new ClippingOp(start, stop));
		final GaeaSamRecord clippedRead = clipRead(clippingOp, runAsserts);
		this.ops = null;
		return clippedRead;
	}

	/**
	 * Soft clip the read to the variable region (from refStart to refStop)
	 * processing also the clipped bases
	 *
	 * @param read
	 *            the read to be clipped
	 * @param refStart
	 *            the beginning of the variant region (inclusive)
	 * @param refStop
	 *            the end of the variant region (inclusive)
	 * @return the read soft clipped to the variant region
	 */
	public static GaeaSamRecord softClipToRegionIncludingClippedBases(final GaeaSamRecord read, final int refStart,
			final int refStop) {
		final int start = read.getUnclippedStart();
		final int stop = start
				+ CigarUtils.countRefBasesBasedOnCigar(read, 0, read.getCigar().getCigarElements().size()) - 1;

		if (start <= refStop && stop >= refStart) {
			if (start < refStart && stop > refStop) {
				return (new ReadClipper(read)).softClipBothEndsByReferenceCoordinates(refStart - 1, refStop + 1);
			} else if (start < refStart) {
				return (new ReadClipper(read)).softClipByReferenceCoordinates(-1, refStart - 1);
			} else if (stop > refStop) {
				return (new ReadClipper(read)).softClipByReferenceCoordinates(refStop + 1, -1);
			}
			return read;
		} else {
			return ReadUtils.emptyRead(read);
		}
	}

	/**
	 * Soft clips both tails of a read. Left tail goes from the beginning to the
	 * 'left' coordinate (inclusive) Right tail goes from the 'right' coordinate
	 * (inclusive) until the end of the read
	 *
	 * @param left
	 *            the coordinate of the last base to be clipped in the left tail
	 *            (inclusive)
	 * @param right
	 *            the coordinate of the first base to be clipped in the right
	 *            tail (inclusive)
	 * @return a new read, without the clipped bases
	 */
	private GaeaSamRecord softClipBothEndsByReferenceCoordinates(final int left, final int right) {
		if (read.isEmpty()) {
			return ReadUtils.emptyRead(read);
		}
		if (left == right) {
			return ReadUtils.emptyRead(read);
		}
		final GaeaSamRecord leftTailRead = softClipByReferenceCoordinates(right, -1);

		// after clipping one tail, it is possible that the consequent hard
		// clipping of adjacent deletions
		// make the left cut index no longer part of the read. In that case,
		// clip the read entirely.
		if (left > leftTailRead.getEnd()) {
			return ReadUtils.emptyRead(read);
		}

		final ReadClipper clipper = new ReadClipper(leftTailRead);
		return clipper.softClipByReferenceCoordinates(-1, left);
	}

	public static GaeaSamRecord softClipBothEndsByReferenceCoordinates(final GaeaSamRecord read, final int left,
			final int right) {
		return (new ReadClipper(read)).softClipBothEndsByReferenceCoordinates(left, right);
	}

	/**
	 * Generic functionality to soft clip a read (is analogous to
	 * hardClipByReferenceCoordinates())
	 *
	 * Note, it REQUIRES you to give the directionality of your soft clip (i.e.
	 * whether you're clipping the left of right tail) by specifying either
	 * refStart < 0 or refStop < 0.
	 *
	 * @param refStart
	 *            first base to clip (inclusive)
	 * @param refStop
	 *            last base to clip (inclusive)
	 * @return a new read, with the soft clipped bases
	 */
	protected GaeaSamRecord softClipByReferenceCoordinates(final int refStart, final int refStop) {
		return clipByReferenceCoordinates(refStart, refStop, ClippingRepresentation.SOFTCLIP_BASES, true);
	}

	/**
	 * Soft clips a read using read coordinates.
	 *
	 * @param start
	 *            the first base to clip (inclusive)
	 * @param stop
	 *            the last base to clip (inclusive)
	 * @return a new read, without the clipped bases
	 */
	private GaeaSamRecord softClipByReadCoordinates(final int start, final int stop) {
		if (read.isEmpty()) {
			return ReadUtils.emptyRead(read);
		}
		if ((start == 0 && stop == read.getReadLength() - 1)) {
			return ReadUtils.emptyRead(read);
		}

		this.addOp(new ClippingOp(start, stop));
		return clipRead(ClippingRepresentation.SOFTCLIP_BASES);
	}

	public static GaeaSamRecord softClipByReadCoordinates(final GaeaSamRecord read, final int start, final int stop) {
		return (new ReadClipper(read)).softClipByReadCoordinates(start, stop);
	}
}
