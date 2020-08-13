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
package org.bgi.flexlab.gaea.data.structure.bam.filter.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.region.Region;

public class MalformedReadFilter implements SamRecordFilter {

	protected boolean checkHasReadGroup(SAMRecord read) {
		if (read.getReadGroup() == null)
			throw new UserException.ReadMissingReadGroup(read);
		return true;
	}

	/**
	 * Check for the case in which the alignment start is inconsistent with the
	 * read unmapped flag.
	 * 
	 * @param read
	 *            The read to validate.
	 * @return true if read start is valid, false otherwise.
	 */
	protected boolean checkInvalidAlignmentStart(SAMRecord read) {
		// read is not flagged as 'unmapped', but alignment start is
		// NO_ALIGNMENT_START
		if (!read.getReadUnmappedFlag()
				&& read.getAlignmentStart() == GaeaSamRecord.NO_ALIGNMENT_START)
			return false;
		// Read is not flagged as 'unmapped', but alignment start is -1
		if (!read.getReadUnmappedFlag() && read.getAlignmentStart() == -1)
			return false;
		return true;
	}

	/**
	 * Check for invalid end of alignments.
	 * 
	 * @param read
	 *            The read to validate.
	 * @return true if read end is valid, false otherwise.
	 */
	protected boolean checkInvalidAlignmentEnd(SAMRecord read) {
		// Alignment aligns to negative number of bases in the reference.
		if (!read.getReadUnmappedFlag() && read.getAlignmentEnd() != -1
				&& (read.getAlignmentEnd() - read.getAlignmentStart() + 1) < 0)
			return false;
		return true;
	}

	/**
	 * Check to ensure that the alignment makes sense based on the contents of
	 * the header.
	 * 
	 * @param header
	 *            The SAM file header.
	 * @param read
	 *            The read to verify.
	 * @return true if alignment agrees with header, false othrewise.
	 */
	protected boolean checkAlignmentDisagreesWithHeader(
			SAMFileHeader header, SAMRecord read) {
		// Read is aligned to nonexistent contig
		if (read.getReferenceIndex() == GaeaSamRecord.NO_ALIGNMENT_REFERENCE_INDEX
				&& read.getAlignmentStart() != GaeaSamRecord.NO_ALIGNMENT_START)
			return false;
		SAMSequenceRecord contigHeader = header.getSequence(read
				.getReferenceIndex());
		// Read is aligned to a point after the end of the contig
		if (!read.getReadUnmappedFlag()
				&& read.getAlignmentStart() > contigHeader.getSequenceLength())
			return false;
		return true;
	}

	/**
	 * Check for inconsistencies between the cigar string and the
	 * 
	 * @param read
	 *            The read to validate.
	 * @return true if cigar agrees with alignment, false otherwise.
	 */
	protected boolean checkCigarDisagreesWithAlignment(SAMRecord read) {
		// Read has a valid alignment start, but the CIGAR string is empty
		if (!read.getReadUnmappedFlag() && read.getAlignmentStart() != -1
				&& read.getAlignmentStart() != GaeaSamRecord.NO_ALIGNMENT_START
				&& read.getAlignmentBlocks().size() < 0)
			return false;
		return true;
	}

	/**
	 * Check if the read has the same number of bases and base qualities
	 * 
	 * @param read
	 *            the read to validate
	 * @return true if they have the same number. False otherwise.
	 */
	protected boolean checkMismatchingBasesAndQuals(SAMRecord read,
			boolean filterMismatchingBaseAndQuals) {
		boolean result;
		if (read.getReadLength() == read.getBaseQualities().length)
			result = true;
		else if (filterMismatchingBaseAndQuals)
			result = false;
		else
			throw new UserException.MalformedBAM(
					read,
					String.format(
							"BAM file has a read with mismatching number of bases and base qualities. Offender: %s [%d bases] [%d quals]",
							read.getReadName(), read.getReadLength(),
							read.getBaseQualities().length));

		return result;
	}

	@Override
	public boolean filter(SAMRecord sam, Region region) {
		return !checkInvalidAlignmentStart(sam)
				|| !checkInvalidAlignmentEnd(sam)
				|| !checkAlignmentDisagreesWithHeader(sam.getHeader(), sam)
				|| !checkHasReadGroup(sam)
				|| !checkMismatchingBasesAndQuals(sam,
						false)
				|| !checkCigarDisagreesWithAlignment(sam);
	}
}
