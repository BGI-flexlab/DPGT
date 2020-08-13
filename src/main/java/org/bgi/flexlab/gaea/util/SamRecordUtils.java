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
package org.bgi.flexlab.gaea.util;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import org.bgi.flexlab.gaea.data.exception.UserException;

import java.util.Arrays;

public class SamRecordUtils {

	public static final String REDUCED_READ_CONSENSUS_TAG = "RR";
	public static final String REDUCED_READ_ORIGINAL_ALIGNMENT_START_SHIFT = "OP";
	public static final String REDUCED_READ_ORIGINAL_ALIGNMENT_END_SHIFT = "OE";

	public static final String BQSR_BASE_INSERTION_QUALITIES = "BI";
	public static final String BQSR_BASE_DELETION_QUALITIES = "BD";

	private final static int UNINITIALIZED = -1;

	public static byte[] getBaseQualities(SAMRecord sam,
			final EventType errorModel) {
		switch (errorModel) {
		case SNP:
			return sam.getBaseQualities();
		case Insertion:
			return getBaseInsertionQualities(sam);
		case Deletion:
			return getBaseDeletionQualities(sam);
		default:
			throw new UserException("Unrecognized Base Recalibration type: "
					+ errorModel);
		}
	}

	public static byte[] getBaseInsertionQualities(SAMRecord sam) {
		byte[] quals = getExistingBaseInsertionQualities(sam);
		if (quals == null) {
			quals = new byte[sam.getBaseQualities().length];
			Arrays.fill(quals, (byte) 45); // 45 is default quality
		}
		return quals;
	}

	public static byte[] getBaseDeletionQualities(SAMRecord sam) {
		byte[] quals = getExistingBaseDeletionQualities(sam);
		if (quals == null) {
			quals = new byte[sam.getBaseQualities().length];
			Arrays.fill(quals, (byte) 45); // the original quality is a flat Q45
		}
		return quals;
	}

	public static byte[] getExistingBaseDeletionQualities(SAMRecord sam) {
		return SAMUtils.fastqToPhred(sam
				.getStringAttribute(BQSR_BASE_DELETION_QUALITIES));
	}

	public static byte[] getExistingBaseInsertionQualities(SAMRecord sam) {
		return SAMUtils.fastqToPhred(sam
				.getStringAttribute(BQSR_BASE_INSERTION_QUALITIES));
	}

	public static int getSoftStart(SAMRecord sam) {
		int softStart = UNINITIALIZED;
		if (softStart == UNINITIALIZED) {
			softStart = sam.getAlignmentStart();
			for (final CigarElement cig : sam.getCigar().getCigarElements()) {
				final CigarOperator op = cig.getOperator();

				if (op == CigarOperator.SOFT_CLIP)
					softStart -= cig.getLength();
				else if (op != CigarOperator.HARD_CLIP)
					break;
			}
		}
		return softStart;
	}

	public static boolean isReducedRead(SAMRecord sam) {
		return getReducedReadCounts(sam) != null;
	}

	public static byte[] getReducedReadCounts(SAMRecord sam) {
		byte[] reducedReadCounts = sam
				.getByteArrayAttribute(REDUCED_READ_CONSENSUS_TAG);

		return reducedReadCounts;
	}

	public static byte getReducedCount(SAMRecord sam, final int index) {
		byte firstCount = getReducedReadCounts(sam)[0];
		byte offsetCount = getReducedReadCounts(sam)[index];
		return (index == 0) ? firstCount : (byte) Math.min(firstCount
				+ offsetCount, Byte.MAX_VALUE);
	}
	
	public static String getReadGroup(SAMRecord record)
	{
		String sam = record.getSAMString();
		String temp;
		int index = sam.indexOf("RG:Z:");
		if(index > 0)
		{
			temp = sam.substring(index);
			index = temp.indexOf("\t");
			if(index < 0) {
				temp = temp.trim();
				index = temp.length();
			}
			temp = temp.substring(0, index);
			String[] splite = temp.split(":");
			return splite[2];
		}else
			return null;
	}
	
	public static boolean isUnmapped(SAMRecord sam) {
		if (sam.getReadUnmappedFlag() || sam.getReferenceIndex() == -1
				|| sam.getAlignmentStart() < 0)
			return true;
		return false;
	}
}
