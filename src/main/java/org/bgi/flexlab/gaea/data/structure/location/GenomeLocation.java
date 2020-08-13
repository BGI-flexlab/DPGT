/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.structure.location;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.input.bed.RegionHdfsParser;
import org.bgi.flexlab.gaea.util.Utils;
import org.bgi.flexlab.gaea.util.Window;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;

import java.io.Serializable;
import java.util.*;

public class GenomeLocation implements Comparable<GenomeLocation>, Comparator<GenomeLocation>, Serializable, Locatable {
	private static final long serialVersionUID = -1554058939692084601L;

	protected final int contigIndex;
	protected final int contigStart;
	protected final int contigEnd;
	protected final String contigName;

	public static final GenomeLocation UNMAPPED = new GenomeLocation((String) null);
	public static final GenomeLocation WHOLE_GENOME = new GenomeLocation("all");

	public static final boolean isUnmapped(GenomeLocation location) {
		return location == UNMAPPED;
	}

	public GenomeLocation(final String contig, final int contigIndex, final int start, final int stop) {
		this.contigName = contig;
		this.contigIndex = contigIndex;
		this.contigStart = start;
		this.contigEnd = stop;
	}
	
	public GenomeLocation(final String contig, final int start, final int stop) {
		this.contigName = contig;
		this.contigIndex = -1;
		this.contigStart = start;
		this.contigEnd = stop;
	}

	public GenomeLocation(final Locatable locatable) {
		this(Utils.nonNull(locatable).getContig(), -1, locatable.getStart(), locatable.getEnd());
	}
	
	public GenomeLocation(final GenomeLocation location) {
		this(Utils.nonNull(location).getContig(), location.getContigIndex(), location.getStart(), location.getEnd());
	}

	/** Unsafe constructor for special constant genome locs */
	private GenomeLocation(final String contig) {
		this(contig, -1, 0, 0);
	}

	public final GenomeLocation getLocation() {
		return this;
	}

	public final GenomeLocation createGenomeLocation(int position) {
		return new GenomeLocation(getContig(), getContigIndex(), position, position);
	}

	public final GenomeLocation createGenomeLocation(int start, int end) {
		return new GenomeLocation(getContig(), getContigIndex(), start, end);
	}
	
	public static final GenomeLocation createGenomeLocation(String contig,int start, int end,int contigLength) {
		return new GenomeLocation(contig, -1, start, Math.min(end, contigLength));
	}
	
	public static final GenomeLocation createGenomeLocation(String contig,int start, int end,int contigStart,int contigLength) {
		return new GenomeLocation(contig, -1, Math.max(start,contigStart), Math.min(end, contigLength));
	}

	public final GenomeLocation getStartLocation() {
		return createGenomeLocation(getStart());
	}

	public final GenomeLocation getStopLocation() {
		return createGenomeLocation(getStop());
	}

	public final String getContig() {
		return this.contigName;
	}

	public final int getContigIndex() {
		return this.contigIndex;
	}

	public final int getStart() {
		return this.contigStart;
	}

	public final int getStop() {
		return getEnd();
	}

	public final String toString() {
		if (GenomeLocation.isUnmapped(this))
			return "unmapped";
		if (throughEndOfContig() && atBeginningOfContig())
			return getContig();
		else if (throughEndOfContig() || getStart() == getStop())
			return String.format("%s:%d", getContig(), getStart());
		else
			return String.format("%s:%d-%d", getContig(), getStart(), getStop());
	}

	private boolean throughEndOfContig() {
		return this.contigEnd == Integer.MAX_VALUE;
	}

	private boolean atBeginningOfContig() {
		return this.contigStart == 1;
	}

	public final boolean disjoint(Locatable that) {
		return contigName.equals(that.getContig()) || this.contigStart > that.getEnd()
				|| that.getStart() > this.contigEnd;
	}

	public final boolean discontinuous(Locatable that) {
		return contigName.equals(that.getContig()) || (this.contigStart - 1) > that.getEnd()
				|| (that.getStart() - 1) > this.contigEnd;
	}

	public final boolean overlaps(Locatable that) {
		return overlapsWithMargin(that,0);
	}
	
	public boolean overlapsWithMargin(final Locatable other, final int margin) {
        if ( margin < 0 ) {
            throw new IllegalArgumentException("given margin is negative: " + margin +
                    "\tfor this: " + toString() + "\tand that: " + (other == null ? "other is null" : other.toString()));
        }
        if ( other == null || other.getContig() == null ) {
            return false;
        }

        return this.contigName.equals(other.getContig()) && this.contigStart <= other.getEnd() + margin && other.getStart() - margin <= this.contigEnd;
    }

	public final boolean contiguous(Locatable that) {
		Utils.nonNull(that);
	    return this.getContig().equals(that.getContig()) && this.getStart() <= that.getEnd() + 1 && that.getStart() <= this.getEnd() + 1;
	}

	/**
	 * Returns a new GenomeLoc that represents the entire span of this and that.
	 * Requires that this and that GenomeLoc are contiguous and both mapped
	 */

	public GenomeLocation merge(GenomeLocation that) throws UserException {
		if (GenomeLocation.isUnmapped(this) || GenomeLocation.isUnmapped(that)) {
			if (!GenomeLocation.isUnmapped(this) || !GenomeLocation.isUnmapped(that))
				throw new UserException("Tried to merge a mapped and an unmapped genome location");
			return UNMAPPED;
		}

		if (!(this.contiguous(that))) {
			throw new UserException("The two genome loc's need to be contiguous");
		}

		return createGenomeLocation(Math.min(getStart(), that.getStart()), Math.max(getStop(), that.getStop()));
	}

	/**
	 * Returns a new GenomeLoc that represents the region between the endpoints
	 * of this and that. Requires that this and that GenomeLoc are both mapped.
	 */
	public GenomeLocation endpointSpan(GenomeLocation that) throws UserException {
		if (GenomeLocation.isUnmapped(this) || GenomeLocation.isUnmapped(that)) {
			throw new UserException("Cannot get endpoint span for unmerged genome locs");
		}

		if (!this.getContig().equals(that.getContig())) {
			throw new UserException("Cannot get endpoint span for genome locs on different contigs");
		}

		return createGenomeLocation(Math.min(getStart(), that.getStart()), Math.max(getStop(), that.getStop()));
	}

	/**
	 * Splits the contig into to regions: [start,split point) and [split point,
	 * end].
	 */
	public GenomeLocation[] split(final int splitPoint) {
		if (splitPoint < getStart() || splitPoint > getStop())
			throw new UserException(String.format(
					"Unable to split contig %s at split point %d; split point is not contained in region.", this,
					splitPoint));
		return new GenomeLocation[] { createGenomeLocation(getStart(), splitPoint - 1),
				createGenomeLocation(splitPoint, getStop()) };
	}

	public GenomeLocation intersect(GenomeLocation that) throws UserException {
		if (GenomeLocation.isUnmapped(this) || GenomeLocation.isUnmapped(that)) {
			if (!GenomeLocation.isUnmapped(this) || !GenomeLocation.isUnmapped(that))
				throw new UserException("Tried to intersect a mapped and an unmapped genome loc");
			return UNMAPPED;
		}

		if (!(this.overlaps(that))) {
			throw new UserException("GenomeLoc::intersect(): The two genome loc's need to overlap");
		}

		return createGenomeLocation(Math.max(getStart(), that.getStart()), Math.min(getStop(), that.getStop()));
	}

	public final List<GenomeLocation> subtract(final GenomeLocation that) {
		if (GenomeLocation.isUnmapped(this) || GenomeLocation.isUnmapped(that)) {
			if (!GenomeLocation.isUnmapped(this) || !GenomeLocation.isUnmapped(that))
				throw new UserException("Tried to intersect a mapped and an unmapped genome loc");
			return Arrays.asList(UNMAPPED);
		}

		if (!(this.overlaps(that))) {
			throw new UserException("GenomeLoc::minus(): The two genome location's need to overlap");
		}

		if (equals(that) || that.contains(this)) {
			return Collections.emptyList();
		} else if (contains(that)) {
			List<GenomeLocation> l = new ArrayList<GenomeLocation>(2);

			int afterStop = this.getStop(), afterStart = that.getStop() + 1;
			int beforeStop = that.getStart() - 1, beforeStart = this.getStart();
			if (afterStop - afterStart >= 0) {
				GenomeLocation after = createGenomeLocation(afterStart, afterStop);
				l.add(after);
			}
			if (beforeStop - beforeStart >= 0) {
				GenomeLocation before = createGenomeLocation(beforeStart, beforeStop);
				l.add(before);
			}
			return l;
		} else {
			GenomeLocation n;
			if (that.getStart() < this.getStart()) {
				n = createGenomeLocation(that.getStop() + 1, this.getStop());
			} else {
				n = createGenomeLocation(this.getStart(), that.getStart() - 1);
			}

			// replace g with the new region
			return Arrays.asList(n);
		}
	}

	public final boolean onSameContig(GenomeLocation that) {
		return (this.contigIndex == that.contigIndex);
	}

	public final int distance(final GenomeLocation that) {
		if (this.onSameContig(that))
			return Math.abs(this.getStart() - that.getStart());
		else
			return Integer.MAX_VALUE;
	}

	public final boolean isBetween(final GenomeLocation left, final GenomeLocation right) {
		return this.compareTo(left) > -1 && this.compareTo(right) < 1;
	}

	/**
	 * Tests whether this contig is completely before contig 'that'.
	 * 
	 */
	public final boolean isBefore(GenomeLocation that) {
		int comparison = this.compareContigs(that);
		return (comparison == -1 || (comparison == 0 && this.getStop() < that.getStart()));
	}

	/**
	 * Tests whether any portion of this contig is before that contig.
	 */
	public final boolean startsBefore(final GenomeLocation that) {
		int comparison = this.compareContigs(that);
		return (comparison == -1 || (comparison == 0 && this.getStart() < that.getStart()));
	}

	/**
	 * Tests whether this contig is completely after contig 'that'.
	 */
	public final boolean isPast(GenomeLocation that) {
		int comparison = this.compareContigs(that);
		return (comparison == 1 || (comparison == 0 && this.getStart() > that.getStop()));
	}

	/**
	 * Return the minimum distance between any pair of bases in this and that
	 * GenomeLocs:
	 */
	public final int minDistance(final GenomeLocation that) {
		if (!this.onSameContig(that))
			return Integer.MAX_VALUE;

		int minDistance;
		if (this.isBefore(that))
			minDistance = distanceFirstStopToSecondStart(this, that);
		else if (that.isBefore(this))
			minDistance = distanceFirstStopToSecondStart(that, this);
		else
			minDistance = 0;

		return minDistance;
	}

	private static int distanceFirstStopToSecondStart(GenomeLocation locFirst, GenomeLocation locSecond) {
		return locSecond.getStart() - locFirst.getStop();
	}

	/**
	 * Check to see whether two genomeLocs are equal. Note that this
	 * implementation ignores the contigInfo object.
	 */
	@Override
	public boolean equals(Object other) {
		if (other == null)
			return false;
		if (other instanceof GenomeLocation) {
			GenomeLocation otherGenomeLoc = (GenomeLocation) other;
			return this.contigIndex == otherGenomeLoc.contigIndex && this.contigStart == otherGenomeLoc.contigStart
					&& this.contigEnd == otherGenomeLoc.contigEnd;
		}
		return false;
	}

	@Override
	public int hashCode() {
		return contigStart << 16 | contigEnd << 4 | contigIndex;
	}

	/**
	 * conpare this genomeLoc's contig to another genome loc
	 */
	public final int compareContigs(GenomeLocation that) {
		if (this.contigIndex == that.contigIndex)
			return 0;
		else if (this.contigIndex > that.contigIndex)
			return 1;
		return -1;
	}

	public boolean endsAt(GenomeLocation that) {
		return (this.compareContigs(that) == 0) && (this.getStop() == that.getStop());
	}

	/**
	 * How many BPs are covered by this locus?
	 */
	public int size() {
		return contigEnd - contigStart + 1;
	}

	/**
	 * reciprocialOverlap: what is the min. percent of gl1 and gl2 covered by
	 * both
	 */
	public final double reciprocialOverlapFraction(final GenomeLocation o) {
		if (overlaps(o))
			return Math.min(overlapPercent(this, o), overlapPercent(o, this));
		else
			return 0.0;
	}

	private final static double overlapPercent(final GenomeLocation gl1, final GenomeLocation gl2) {
		return (1.0 * gl1.intersect(gl2).size()) / gl1.size();
	}

	public long sizeOfOverlap(final GenomeLocation that) {
		return (this.overlaps(that) ? Math.min(getStop(), that.getStop()) - Math.max(getStart(), that.getStart()) + 1L
				: 0L);
	}

	@Override
	public int compare(GenomeLocation g1, GenomeLocation g2) {
		int cmp = g1.getContigIndex() - g2.getContigIndex();
		if (cmp == 0) {
			cmp = g1.getStart() - g2.getStart();
			if (cmp == 0) {
				cmp = g1.getStop() - g2.getStop();
			}
		}
		return cmp;
	}

	@Override
	public int compareTo(GenomeLocation that) {
		int result = 0;
		if (this == that) {
			result = 0;
		} else if (GenomeLocation.isUnmapped(this)) {
			result = 1;
		} else if (GenomeLocation.isUnmapped(that)) {
			result = -1;
		} else {
			result = compare(this, that);
		}
		return result;
	}

	public static final int compareLocatables(GenomeLocation first, GenomeLocation second,
			SAMSequenceDictionary dictionary) {
		Utils.nonNull(first);
		Utils.nonNull(second);
		Utils.nonNull(dictionary);

		int result = 0;
		if (first != second) {
			// compare the contigs
			result = compareContigs(first, second, dictionary);
			if (result == 0) {
				// compare start position
				result = Integer.compare(first.getStart(), second.getStart());
				if (result == 0) {
					// compare end position
					result = Integer.compare(first.getStop(), second.getStop());
				}
			}
		}
		return result;
	}

	public static int compareContigs(final GenomeLocation first, final GenomeLocation second,
			final SAMSequenceDictionary dictionary) {
		Utils.nonNull(first);
		Utils.nonNull(second);
		Utils.nonNull(dictionary);

		final int firstRefIndex = dictionary.getSequenceIndex(first.getContig());
		final int secondRefIndex = dictionary.getSequenceIndex(second.getContig());
		if (firstRefIndex == -1 || secondRefIndex == -1) {
			throw new IllegalArgumentException(
					"Can't do comparison because Locatables' contigs not found in sequence dictionary");
		}
		// compare the contigs
		return Integer.compare(firstRefIndex, secondRefIndex);
	}

	@Override
	public int getEnd() {
		return this.contigEnd;
	}

	public long getGA4GHEnd() {
		return this.contigEnd;
	}

	public long getGA4GHStart() {
		return contigStart - 1;
	}

	public boolean contains(final Locatable other) {
		if (other == null || other.getContig() == null) {
			return false;
		}

		return this.contigName.equals(other.getContig()) && this.contigStart <= other.getStart()
				&& this.contigEnd >= other.getEnd();
	}
	
	public GenomeLocation expandWithinContig(final int padding, final int contigLength) {
		Utils.validateArg(padding >= 0, "padding must be >= 0");
		final int boundedStart = Math.max(1, contigStart - padding);
		final int boundedStop = Math.min(contigLength, contigEnd + padding);
		return createGenomeLocation(boundedStart, boundedStop);
	}
	
	public GenomeLocation expandWithinContig(final int padding, final SAMSequenceDictionary sequenceDictionary) {
		Utils.nonNull(sequenceDictionary);
		final SAMSequenceRecord contigRecord = sequenceDictionary.getSequence(contigName);
		return expandWithinContig(padding, contigRecord.getSequenceLength());
	}
	
	public GenomeLocation spanWith(final Locatable other) {
		Utils.nonNull(other);
		Utils.validateArg(this.getContig().equals(other.getContig()),
				"Cannot get span for intervals on different contigs");
		return new GenomeLocation(contigName, Math.min(contigStart, other.getStart()), Math.max(contigEnd, other.getEnd()));
	}
	
	public GenomeLocation mergeWithContiguous(final Locatable that) {
		Utils.nonNull(that);
		if (!this.contiguous(that)) {
			throw new UserException("The two intervals need to be contiguous: " + this + " " + that);
		}

		return new GenomeLocation(getContig(), Math.min(getStart(), that.getStart()),
				Math.max(getEnd(), that.getEnd()));
	}
	
	public static void getGenomeLocationFromWindow(List<GenomeLocation> intervals,Window win,RegionHdfsParser region){
		if (!intervals.isEmpty())
			intervals.clear();

		String chr = win.getContigName();
		int start = win.getStart();
		int index = win.getChrIndex();
		int stop = win.getStop();
		if (region == null) {
			intervals.add(new GenomeLocation(chr, index, start, stop));
			return ;
		}
		int length = stop - start;

		BitSet bitSet = new BitSet(length);
		for (int i = start; i < stop; i++) {
			if (region.isPositionInRegion(chr, i))
				bitSet.set(i - start);
		}

		int intervalStart = -1;

		for (int i = 0; i < length; i++) {
			if (bitSet.get(i)) {
				if (intervalStart == -1) {
					intervalStart = i + start;
				}
			} else {
				if (intervalStart != -1) {
					intervals.add(new GenomeLocation(chr, index, intervalStart, i + start));
					intervalStart = -1;
				}
			}
		}
		bitSet = null;

		if (intervalStart != -1)
			intervals.add(new GenomeLocation(chr, index, intervalStart, stop));
	} 
}
