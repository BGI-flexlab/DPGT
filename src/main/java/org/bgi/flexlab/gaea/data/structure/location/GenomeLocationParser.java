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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.bgi.flexlab.gaea.data.exception.UserException;

public class GenomeLocationParser {
	private SAMSequenceDictionary SINGLE_MASTER_SEQUENCE_DICTIONARY;

	private final ThreadLocal<CachingSequenceDictionary> contigInfoPerThread = new ThreadLocal<CachingSequenceDictionary>();

	private CachingSequenceDictionary getContigInfo() {
		if (contigInfoPerThread.get() == null) {
			// initialize for this thread
			contigInfoPerThread.set(new CachingSequenceDictionary(SINGLE_MASTER_SEQUENCE_DICTIONARY));
		}

		assert contigInfoPerThread.get() != null;

		return contigInfoPerThread.get();
	}

	private final class CachingSequenceDictionary {
		final private SAMSequenceDictionary dict;

		// cache
		SAMSequenceRecord lastSSR = null;
		String lastContig = "";
		int lastIndex = -1;

		// @Requires({"dict != null", "dict.size() > 0"})
		public CachingSequenceDictionary(SAMSequenceDictionary dict) {
			this.dict = dict;
		}

		// @Ensures("result > 0")
		public final int getNSequences() {
			return dict.size();
		}

		// @Requires("contig != null")
		public final boolean hasContig(final String contig) {
			return contig.equals(lastContig) || dict.getSequence(contig) != null;
		}

		// @Requires("index >= 0")
		public final boolean hasContig(final int index) {
			return lastIndex == index || dict.getSequence(index) != null;
		}

		// @Requires("contig != null")
		// @Ensures("result != null")
		public final SAMSequenceRecord getSequence(final String contig) {
			if (isCached(contig))
				return lastSSR;
			else
				return updateCache(contig, -1);
		}

		// @Requires("index >= 0")
		// @Ensures("result != null")
		public final SAMSequenceRecord getSequence(final int index) {
			if (isCached(index))
				return lastSSR;
			else
				return updateCache(null, index);
		}

		// @Requires("contig != null")
		// @Ensures("result >= 0")
		public final int getSequenceIndex(final String contig) {
			if (!isCached(contig)) {
				updateCache(contig, -1);
			}

			return lastIndex;
		}

		// @Requires({"contig != null", "lastContig != null"})
		private boolean isCached(final String contig) {
			return lastContig.equals(contig);
		}

		// @Requires({"lastIndex != -1", "index >= 0"})
		private boolean isCached(final int index) {
			return lastIndex == index;
		}

		/**
		 * The key algorithm. Given a new record, update the last used record,
		 * contig name, and index.
		 */
		private SAMSequenceRecord updateCache(final String contig, int index) {
			SAMSequenceRecord rec = contig == null ? dict.getSequence(index) : dict.getSequence(contig);
			if (rec == null) {
				throw new UserException("BUG: requested unknown contig=" + contig + " index=" + index);
			} else {
				lastSSR = rec;
				lastContig = rec.getSequenceName();
				lastIndex = rec.getSequenceIndex();
				return rec;
			}
		}

	}

	public GenomeLocationParser() {

	}

	public GenomeLocationParser(SAMSequenceDictionary seqDict) {
		if (seqDict == null) { // we couldn't load the reference dictionary
			throw new UserException("Failed to load reference dictionary");
		}

		SINGLE_MASTER_SEQUENCE_DICTIONARY = seqDict;
	}

	/**
	 * Determines whether the given contig is valid with respect to the sequence
	 * dictionary already installed in the GenomeLoc.
	 */
	public final boolean contigIsInDictionary(String contig) {
		return contig != null && getContigInfo().hasContig(contig);
	}

	public final boolean indexIsInDictionary(final int index) {
		return index >= 0 && getContigInfo().hasContig(index);
	}

	/**
	 * get the contig's SAMSequenceRecord
	 */
	public final SAMSequenceRecord getContigInfo(final String contig) {
		if (contig == null || !contigIsInDictionary(contig))
			throw new UserException(String.format(
					"Contig %s given as location, but this contig isn't present in the Fasta sequence dictionary",
					contig));
		return getContigInfo().getSequence(contig);
	}

	public final SAMSequenceRecord getContigInfo(final int index) {
		if (!indexIsInDictionary(index))
			throw new UserException(String.format(
					"Contig index %s given as location, but this contig isn't present in the Fasta sequence dictionary",
					index));
		return getContigInfo().getSequence(index);
	}

	/**
	 * Returns the contig index of a specified string version of the contig
	 */
	public final int getContigIndex(final String contig) {
		return getContigInfo(contig).getSequenceIndex();
	}

	protected int getContigIndexWithoutException(final String contig) {
		if (contig == null || !getContigInfo().hasContig(contig))
			return -1;
		return getContigInfo().getSequenceIndex(contig);
	}

	/**
	 * Return the master sequence dictionary used within this GenomeLocParser
	 */
	public final SAMSequenceDictionary getContigs() {
		return getContigInfo().dict;
	}

	/**
	 * create a genome location, given the contig name, start, and stop
	 */
	public GenomeLocation createGenomeLocation(String contig, final int start, final int stop) {
		return createGenomeLocation(contig, getContigIndex(contig), start, stop);
	}

	public GenomeLocation createGenomeLocation(String contig, final int start, final int stop, boolean mustBeOnReference) {
		return createGenomeLocation(contig, getContigIndex(contig), start, stop, mustBeOnReference);
	}

	public GenomeLocation createGenomeLocation(String contig, int index, final int start, final int stop) {
		return createGenomeLocation(contig, index, start, stop, false);
	}

	public GenomeLocation createGenomeLocation(String contig, int index, final int start, final int stop,
			boolean mustBeOnReference) {
		validateGenomeLocation(contig, index, start, stop, mustBeOnReference, true);
		return new GenomeLocation(contig, index, start, stop);
	}

	/**
	 * validate a position or interval on the genome as valid
	 */
	private boolean validateGenomeLocation(String contig, int contigIndex, int start, int stop,
			boolean mustBeOnReference, boolean exceptOnError) {
		if (!getContigInfo().hasContig(contig))
			return vglHelper(exceptOnError, String.format("Unknown contig %s", contig));

		if (stop < start)
			return vglHelper(exceptOnError,
					String.format("The stop position %d is less than start %d in contig %s", stop, start, contig));

		if (contigIndex < 0)
			return vglHelper(exceptOnError, String.format("The contig index %d is less than 0", contigIndex));

		if (contigIndex >= getContigInfo().getNSequences())
			return vglHelper(exceptOnError,
					String.format("The contig index %d is greater than the stored sequence count (%d)", contigIndex,
							getContigInfo().getNSequences()));

		if (mustBeOnReference) {
			if (start < 1)
				return vglHelper(exceptOnError, String.format("The start position %d is less than 1", start));

			if (stop < 1)
				return vglHelper(exceptOnError, String.format("The stop position %d is less than 1", stop));

			int contigSize = getContigInfo().getSequence(contigIndex).getSequenceLength();
			if (start > contigSize || stop > contigSize)
				return vglHelper(exceptOnError, String.format(
						"The genome loc coordinates %d-%d exceed the contig size (%d)", start, stop, contigSize));
		}

		// we passed
		return true;
	}

	public boolean isValidGenomeLocation(String contig, int start, int stop, boolean mustBeOnReference) {
		return validateGenomeLocation(contig, getContigIndexWithoutException(contig), start, stop, mustBeOnReference,
				false);
	}

	public boolean isValidGenomeLocation(String contig, int start, int stop) {
		return validateGenomeLocation(contig, getContigIndexWithoutException(contig), start, stop, true, false);
	}

	private boolean vglHelper(boolean exceptOnError, String msg) {
		if (exceptOnError)
			throw new UserException("Parameters to GenomeLocParser are incorrect:" + msg);
		else
			return false;
	}

	/**
	 * parse a genome interval, from a location string
	 */
	public GenomeLocation parseGenomeLocation(final String str) {

		String contig = null;
		int start = 1;
		int stop = -1;

		final int colonIndex = str.indexOf(":");
		if (colonIndex == -1) {
			contig = str.substring(0, str.length()); // chr1
			stop = Integer.MAX_VALUE;
		} else {
			contig = str.substring(0, colonIndex);
			final int dashIndex = str.indexOf('-', colonIndex);
			try {
				if (dashIndex == -1) {
					if (str.charAt(str.length() - 1) == '+') {
						start = parsePosition(str.substring(colonIndex + 1, str.length() - 1)); // chr:1+
						stop = Integer.MAX_VALUE;
					} else {
						start = parsePosition(str.substring(colonIndex + 1)); // chr1:1
						stop = start;
					}
				} else {
					start = parsePosition(str.substring(colonIndex + 1, dashIndex)); // chr1:1-1
					stop = parsePosition(str.substring(dashIndex + 1));
				}
			} catch (Exception e) {
				throw new UserException("Failed to parse Genome Location string: " + str, e);
			}
		}

		// is the contig valid?
		if (!contigIsInDictionary(contig))
			throw new UserException("Contig '" + contig
					+ "' does not match any contig in the GATK sequence dictionary derived from the reference; "
					+ "are you sure you are using the correct reference fasta file?");

		if (stop == Integer.MAX_VALUE)
			// lookup the actually stop position!
			stop = getContigInfo(contig).getSequenceLength();

		return createGenomeLocation(contig, getContigIndex(contig), start, stop, true);
	}

	/**
	 * Parses a number like 1,000,000 into a long.
	 */
	private int parsePosition(final String pos) {
		if (pos.indexOf('-') != -1) {
			throw new NumberFormatException("Position: '" + pos + "' can't contain '-'.");
		}

		if (pos.indexOf(',') != -1) {
			final StringBuilder buffer = new StringBuilder();
			for (int i = 0; i < pos.length(); i++) {
				final char c = pos.charAt(i);

				if (c == ',') {
					continue;
				} else if (c < '0' || c > '9') {
					throw new NumberFormatException("Position: '" + pos + "' contains invalid chars.");
				} else {
					buffer.append(c);
				}
			}
			return Integer.parseInt(buffer.toString());
		} else {
			return Integer.parseInt(pos);
		}
	}

	/**
	 * create a genome loc, given a read. If the read is unmapped, *and* yet the
	 * read has a contig and start position, then a GenomeLoc is returned for
	 * contig:start-start, otherwise and UNMAPPED GenomeLoc is returned.
	 */
	public GenomeLocation createGenomeLocation(final SAMRecord read) {
		if (read.getReadUnmappedFlag() && read.getReferenceIndex() == -1)
			// read is unmapped and not placed anywhere on the genome
			return GenomeLocation.UNMAPPED;
		else {
			// Use Math.max to ensure that end >= start (Picard assigns the end
			// to reads that are entirely within an insertion as start-1)
			int end = read.getReadUnmappedFlag() ? read.getAlignmentStart()
					: Math.max(read.getAlignmentEnd(), read.getAlignmentStart());
			return createGenomeLocation(read.getReferenceName(), read.getReferenceIndex(), read.getAlignmentStart(),
					end, false);
		}
	}

	public GenomeLocation createGenomeLocation(final Feature feature) {
		return createGenomeLocation(feature.getContig(), feature.getStart(), feature.getEnd());
	}

	/**
	 * Creates a GenomeLoc corresponding to the variant context vc. If
	 * includeSymbolicEndIfPossible is true, and VC is a symbolic allele the end
	 * of the created genome loc will be the value of the END info field key, if
	 * it exists, or vc.getEnd() if not.
	 */
	public GenomeLocation createGenomeLocation(final VariantContext vc, boolean includeSymbolicEndIfPossible) {
		if (includeSymbolicEndIfPossible && vc.isSymbolic()) {
			int end = vc.getAttributeAsInt(VCFConstants.END_KEY, vc.getEnd());
			return createGenomeLocation(vc.getContig(), vc.getStart(), end);
		} else
			return createGenomeLocation(vc.getContig(), vc.getStart(), vc.getEnd());
	}

	public GenomeLocation createGenomeLocation(final VariantContext vc) {
		return createGenomeLocation(vc, false);
	}

	/**
	 * create a new genome loc, given the contig name, and a single position.
	 * Must be on the reference
	 */
	public GenomeLocation createGenomeLocation(final String contig, final int pos) {
		return createGenomeLocation(contig, getContigIndex(contig), pos, pos);
	}

	public GenomeLocation createGenomeLocation(final int index, final int pos) {
		return createGenomeLocation(getContigInfo(index).getSequenceName(), index, pos, pos);
	}

	/**
	 * create a new genome loc from an existing loc, with a new start position
	 */
	public GenomeLocation setStart(GenomeLocation loc, int start) {
		return createGenomeLocation(loc.getContig(), loc.getContigIndex(), start, loc.getStop());
	}

	public GenomeLocation setStop(GenomeLocation loc, int stop) {
		return createGenomeLocation(loc.getContig(), loc.getContigIndex(), loc.getStart(), stop);
	}

	/**
	 * return a new genome location, with an incremented position
	 */
	public GenomeLocation incPos(GenomeLocation loc) {
		return incPos(loc, 1);
	}

	/**
	 * return a new genome location, with an incremented position
	 */
	public GenomeLocation incPos(GenomeLocation loc, int by) {
		return createGenomeLocation(loc.getContig(), loc.getContigIndex(), loc.getStart() + by, loc.getStop() + by);
	}

	/**
	 * Creates a GenomeLocation than spans the entire contig.
	 */
	public GenomeLocation createOverEntireContig(String contigName) {
		SAMSequenceRecord contig = getContigInfo().getSequence(contigName);
		return createGenomeLocation(contigName, contig.getSequenceIndex(), 1, contig.getSequenceLength(), true);
	}

	/**
	 * Creates a location to the left (starting at the location start + 1) of
	 * maxBasePairs size.
	 */
	public GenomeLocation createGenomeLocAtStart(GenomeLocation loc, int maxBasePairs) {
		if (GenomeLocation.isUnmapped(loc))
			return null;
		String contigName = loc.getContig();
		SAMSequenceRecord contig = getContigInfo().getSequence(contigName);
		int contigIndex = contig.getSequenceIndex();

		int start = loc.getStart() - maxBasePairs;
		int stop = loc.getStart() - 1;

		if (start < 1)
			start = 1;
		if (stop < 1)
			return null;

		return createGenomeLocation(contigName, contigIndex, start, stop, true);
	}

	/**
	 * Creates a location padded in both directions by maxBasePairs size (if
	 * possible).
	 */
	public GenomeLocation createPaddedGenomeLoc(final GenomeLocation loc, final int padding) {
		if (GenomeLocation.isUnmapped(loc))
			return loc;
		final String contigName = loc.getContig();
		final SAMSequenceRecord contig = getContigInfo().getSequence(contigName);
		final int contigIndex = contig.getSequenceIndex();
		final int contigLength = contig.getSequenceLength();

		final int start = Math.max(1, loc.getStart() - padding);
		final int stop = Math.min(contigLength, loc.getStop() + padding);

		return createGenomeLocation(contigName, contigIndex, start, stop, true);
	}

	/**
	 * Creates a location to the right (starting at the loc stop + 1) of
	 * maxBasePairs size.
	 */
	public GenomeLocation createGenomeLocAtStop(GenomeLocation loc, int maxBasePairs) {
		if (GenomeLocation.isUnmapped(loc))
			return null;
		String contigName = loc.getContig();
		SAMSequenceRecord contig = getContigInfo().getSequence(contigName);
		int contigIndex = contig.getSequenceIndex();
		int contigLength = contig.getSequenceLength();

		int start = loc.getStop() + 1;
		int stop = loc.getStop() + maxBasePairs;

		if (start > contigLength)
			return null;
		if (stop > contigLength)
			stop = contigLength;

		return createGenomeLocation(contigName, contigIndex, start, stop, true);
	}
}
