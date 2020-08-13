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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import java.util.Arrays;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaCigar;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;

public class AlignmentUtil {

	public static int mismatchQualityCount(GaeaSamRecord read, byte[] ref, int posOnRef) {
		return mismatchQualityCount(read, ref, posOnRef, 0, read.getReadLength());
	}

	public static int mismatchQualityCount(GaeaSamRecord read, byte[] ref, int posOnRef, int posOnRead,
			int baseLength) {
		int mismatchQualitySum = 0;

		int readIndex = 0;
		int endIndex = posOnRead + baseLength - 1;

		byte[] readSeq = read.getReadBases();
		byte[] qualities = read.getBaseQualities();

		for (CigarElement element : read.getCigar().getCigarElements()) {
			int length = element.getLength();

			switch (element.getOperator()) {
			case M:
				for (int i = 0; i < length; i++, readIndex++, posOnRef++) {
					if (posOnRef >= ref.length || readIndex > endIndex) {
						return mismatchQualitySum;
					}
					if (readIndex < posOnRead)
						continue;
					byte readBase = readSeq[readIndex];
					byte refBase = ref[posOnRef];

					if (readBase != refBase)
						mismatchQualitySum += qualities[readIndex];
				}
				break;
			case X:
				for (int i = 0; i < length; i++, readIndex++, posOnRef++) {
					if (posOnRef > ref.length || readIndex > endIndex) {
						return mismatchQualitySum;
					}
					if (readIndex < posOnRead)
						continue;

					mismatchQualitySum += qualities[readIndex];
				}
				break;
			case EQ:
				readIndex += length;
				posOnRef += length;
				break;
			case I:
			case S:
				readIndex += length;
				break;
			case D:
			case N:
				posOnRef += length;
				break;
			default:
				break;
			}
		}

		return mismatchQualitySum;
	}

	public static byte[] createStringByIndel(Cigar cigar, int indexOfIndel, byte[] ref, byte[] read, int refIndex,
			int readIndex) {
		CigarElement element = cigar.getCigarElement(indexOfIndel);
		int indelLength = element.getLength();

		int i;
		int totalRefBaseCount = 0;
		for (i = 0; i < indexOfIndel; i++) {
			CigarElement ce = cigar.getCigarElement(i);
			switch (ce.getOperator()) {
			case X:
			case EQ:
			case M:
				refIndex += ce.getLength();
				readIndex += ce.getLength();
				totalRefBaseCount += ce.getLength();
				break;
			case S:
				readIndex += ce.getLength();
				break;
			case N:
				refIndex += ce.getLength();
				totalRefBaseCount += ce.getLength();
				break;
			default:
				break;
			}
		}

		/* CigarOperator.I needn't be changed */
		if (element.getOperator() == CigarOperator.D && (indelLength + totalRefBaseCount > ref.length)) {
			indelLength -= (indelLength + totalRefBaseCount - ref.length);
		}

		int refLength = ref.length + (indelLength * (element.getOperator() == CigarOperator.D ? -1 : 1));
		byte[] alt = new byte[refLength];

		if (refIndex > alt.length || refIndex > ref.length)
			return null;

		System.arraycopy(ref, 0, alt, 0, refIndex);

		int currPos = refIndex;

		if (element.getOperator() == CigarOperator.D) {
			refIndex += indelLength;
		} else {
			System.arraycopy(read, readIndex, alt, currPos, indelLength);
			currPos += indelLength;
		}

		if (currPos >= alt.length || (ref.length - refIndex > alt.length - currPos))
			return null;

		System.arraycopy(ref, refIndex, alt, currPos, ref.length - refIndex);

		return alt;
	}

	public static Cigar leftAlignIndel(Cigar originCigar, final byte[] refSeq, final byte[] readSeq, final int refIndex,
			final int readIndex) {
		int indexOfIndel = GaeaCigar.firstIndexOfIndel(originCigar);

		if (indexOfIndel < 1)
			return originCigar;

		final int indelLength = originCigar.getCigarElement(indexOfIndel).getLength();

		byte[] alt = createStringByIndel(originCigar, indexOfIndel, refSeq, readSeq, refIndex, readIndex);
		if (alt == null)
			return originCigar;

		Cigar newCigar = originCigar;
		for (int i = 0; i < indelLength; i++) {
			newCigar = GaeaCigar.moveCigarLeft(newCigar, indexOfIndel, 1);

			if (newCigar == null)
				return originCigar;

			byte[] newAlt = createStringByIndel(newCigar, indexOfIndel, refSeq, readSeq, refIndex, readIndex);

			boolean reachedEndOfRead = GaeaCigar.cigarHasZeroSizeElement(newCigar);

			if (Arrays.equals(alt, newAlt)) {
				originCigar = newCigar;
				i = -1;
				if (reachedEndOfRead)
					originCigar = GaeaCigar.cleanCigar(originCigar);
			}

			if (reachedEndOfRead)
				break;
		}

		return originCigar;
	}

	public static int[] referencePositions(Cigar cigar, int start, int readLength) {
		int[] positions = new int[readLength];
		Arrays.fill(positions, 0);

		int readIndex = 0;
		int referenceIndex = start;

		for (CigarElement element : cigar.getCigarElements()) {
			int length = element.getLength();
			CigarOperator op = element.getOperator();

			switch (op) {
			case M:
			case EQ:
			case X:
				for (int i = 0; i < length; i++) {
					positions[readIndex++] = referenceIndex++;
				}
				break;
			case S:
			case I:
				readIndex += length;
				break;
			case D:
			case N:
				referenceIndex += length;
				break;
			default:
				break;
			}
		}

		return positions;
	}

	public static int[] readOffsets(Cigar cigar, int start, int end) {
		int length = end - start + 1;
		int[] positions = new int[length];
		Arrays.fill(positions, 0);

		CigarState state = new CigarState();
		state.parseCigar(cigar.toString());

		for (int i = start; i <= end; i++) {
			positions[i - start] = state.resolveCigar(i, start);
		}

		return positions;
	}
}
