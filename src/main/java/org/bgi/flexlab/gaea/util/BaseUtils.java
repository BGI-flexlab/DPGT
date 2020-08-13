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

import java.util.Arrays;

public class BaseUtils {
	public final static byte A = (byte) 'A';
	public final static byte C = (byte) 'C';
	public final static byte G = (byte) 'G';
	public final static byte T = (byte) 'T';

	public final static byte N = (byte) 'N';
	public final static byte D = (byte) 'D';

	public final static byte[] BASES = { A, C, T, G };
	public final static byte[] EXTENDED_BASES = { A, C, T, G, N, D };

	private static byte[] complementBaseTable = { 2, 3, 0, 1 };

	public enum Base {
		A('A', 0), C('C', 1), T('T', 2), G('G', 3);

		byte b;
		int index;

		private Base(char base, int index) {
			this.b = (byte) base;
			this.index = index;
		}

		public byte getBase() {
			return b;
		}

		public char getBaseAsChar() {
			return (char) b;
		}

		public int getIndex() {
			return index;
		}

		public boolean sameBase(byte o) {
			return b == o;
		}

		public boolean sameBase(char o) {
			return b == (byte) o;
		}

		public boolean sameBase(int i) {
			return index == i;
		}
	}

	static private final int[] baseIndexMap = new int[256];
	static {
		Arrays.fill(baseIndexMap, -1);
		baseIndexMap['A'] = 0;
		baseIndexMap['a'] = 0;
		baseIndexMap['*'] = 0;
		baseIndexMap['C'] = 1;
		baseIndexMap['c'] = 1;
		baseIndexMap['T'] = 2;
		baseIndexMap['t'] = 2;
		baseIndexMap['G'] = 3;
		baseIndexMap['g'] = 3;
	}

	public static final byte DELETION_INDEX = 5;
	public static final byte NO_CALL_INDEX = 4; // (this is 'N')

	public static final int aIndex = BaseUtils
			.simpleBaseToBaseIndex((byte) 'A');
	public static final int cIndex = BaseUtils
			.simpleBaseToBaseIndex((byte) 'C');
	public static final int gIndex = BaseUtils
			.simpleBaseToBaseIndex((byte) 'G');
	public static final int tIndex = BaseUtils
			.simpleBaseToBaseIndex((byte) 'T');

	public enum BaseSubstitutionType {
		TRANSITION, // A <-> G or C <-> T
		TRANSVERSION // others
	}

	/**
	 * Returns the base substitution type of the 2 state SNP
	 */
	public static BaseSubstitutionType SNPSubstitutionType(byte base1,
			byte base2) {
		BaseSubstitutionType t = isTransition(base1, base2) ? BaseSubstitutionType.TRANSITION
				: BaseSubstitutionType.TRANSVERSION;
		return t;
	}

	/**
	 * A<->G or C<->T is transition
	 */
	public static boolean isTransition(byte base1, byte base2) {
		int b1 = simpleBaseToBaseIndex(base1);
		int b2 = simpleBaseToBaseIndex(base2);
		return ((b1 ^ b2) == 2);
	}

	public static boolean isTransversion(byte base1, byte base2) {
		return !isTransition(base1, base2);
	}

	/**
	 * A->T or C->G is complement
	 * 
	 * @param base
	 * @return
	 */
	public static byte getComplementBase(byte base) {
		return complementBaseTable[base];
	}

	public static byte getBinaryBase(byte charBase) {
		byte base;
		if (charBase == 'N' || charBase == 'n') {
			base = 4;
		} else {
			base = (byte) ((charBase >> 1) & 0x03);
		}
		return base;
	}

	static public boolean basesAreEqual(byte base1, byte base2) {
		return simpleBaseToBaseIndex(base1) == simpleBaseToBaseIndex(base2);
	}

	static public boolean extendedBasesAreEqual(byte base1, byte base2) {
		return extendedBaseToBaseIndex(base1) == extendedBaseToBaseIndex(base2);
	}

	static public boolean containsBase(final byte[] bases, final byte base) {
		for (final byte b : bases) {
			if (b == base)
				return true;
		}
		return false;
	}

	/**
	 * Converts a IUPAC nucleotide code to a pair of bases
	 */
	@Deprecated
	static public char[] iupacToBases(char code) {
		char[] bases = new char[2];
		switch (code) {
		case '*': // the wildcard character counts as an A
		case 'A':
		case 'a':
			bases[0] = bases[1] = 'A';
			break;
		case 'C':
		case 'c':
			bases[0] = bases[1] = 'C';
			break;
		case 'G':
		case 'g':
			bases[0] = bases[1] = 'G';
			break;
		case 'T':
		case 't':
			bases[0] = bases[1] = 'T';
			break;
		case 'R':
		case 'r':
			bases[0] = 'A';
			bases[1] = 'G';
			break;
		case 'Y':
		case 'y':
			bases[0] = 'C';
			bases[1] = 'T';
			break;
		case 'S':
		case 's':
			bases[0] = 'G';
			bases[1] = 'C';
			break;
		case 'W':
		case 'w':
			bases[0] = 'A';
			bases[1] = 'T';
			break;
		case 'K':
		case 'k':
			bases[0] = 'G';
			bases[1] = 'T';
			break;
		case 'M':
		case 'm':
			bases[0] = 'A';
			bases[1] = 'C';
			break;
		default:
			bases[0] = bases[1] = 'N';
		}
		return bases;
	}

	/**
	 * Converts a simple base to a base index
	 */
	static public int simpleBaseToBaseIndex(byte base) {
		return baseIndexMap[base];
	}

	static public int extendedBaseToBaseIndex(byte base) {
		switch (base) {
		case 'd':
		case 'D':
			return DELETION_INDEX;
		case 'n':
		case 'N':
			return NO_CALL_INDEX;

		default:
			return simpleBaseToBaseIndex(base);
		}
	}

	static public boolean isRegularBase(final byte base) {
		return simpleBaseToBaseIndex(base) != -1;
	}

	static public boolean isRegularAndNotEqualBase(final byte readBase,
			final byte refBase) {
		if (isRegularBase(readBase) && isRegularBase(refBase)
				&& (readBase != refBase))
			return true;
		return false;
	}

	static public boolean isAllRegularBases(final byte[] bases) {
		for (final byte base : bases) {
			if (!isRegularBase(base)) {
				return false;
			}
		}
		return true;
	}

	static public boolean isNBase(byte base) {
		return base == 'N' || base == 'n';
	}

	/**
	 * Converts a base index to a simple base
	 */
	static public byte baseIndexToSimpleBase(int baseIndex) {
		switch (baseIndex) {
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'T';
		case 3:
			return 'G';
		default:
			return '.';
		}
	}

	/**
	 * Converts a base index to a base index representing its cross-talk partner
	 */
	static public int crossTalkPartnerIndex(int baseIndex) {
		switch (baseIndex) {
		case 0:
			return 1; // A -> C
		case 1:
			return 0; // C -> A
		case 2:
			return 3; // T -> G
		case 3:
			return 2; // G -> T
		default:
			return -1;
		}
	}

	/**
	 * Return the complement of a base index.
	 */
	static public byte complementIndex(int baseIndex) {
		switch (baseIndex) {
		case 0:
			return 2; // a -> t
		case 1:
			return 3; // c -> g
		case 3:
			return 1; // g -> c
		case 2:
			return 0; // t -> a
		default:
			return -1; // wtf?
		}
	}

	/**
	 * Return the complement (A <-> T or C <-> G) of a base, or the specified
	 * base if it can't be complemented (i.e. an ambiguous base).
	 */
	static public byte simpleComplement(byte base) {
		switch (base) {
		case 'A':
		case 'a':
			return 'T';
		case 'C':
		case 'c':
			return 'G';
		case 'G':
		case 'g':
			return 'C';
		case 'T':
		case 't':
			return 'A';
		default:
			return base;
		}
	}

	/**
	 * Reverse complement a byte array of bases (that is, chars casted to bytes,
	 * *not* base indices in byte form)
	 */
	static public byte[] simpleReverseComplement(byte[] bases) {
		int length = bases.length;
		byte[] rcbases = new byte[length];

		for (int i = 0; i < length; i++) {
			rcbases[i] = simpleComplement(bases[length - 1 - i]);
		}

		return rcbases;
	}

	/**
	 * Complement a byte array of bases (that is, chars casted to bytes, *not*
	 * base indices in byte form)
	 */
	static public byte[] simpleComplement(byte[] bases) {
		int length = bases.length;
		byte[] rcbases = new byte[length];

		for (int i = 0; i < length; i++) {
			rcbases[i] = simpleComplement(bases[i]);
		}

		return rcbases;
	}

	/**
	 * Returns the uppercased version of the bases
	 */
	static public byte[] convertToUpperCase(final byte[] bases) {
		for (int i = 0; i < bases.length; i++) {
			if ((char) bases[i] >= 'a') {
				bases[i] = toUpperCaseBase(bases[i]);
			}
		}
		return bases;
	}

	static public byte toUpperCaseBase(final byte base) {
		switch (base) {
		case 'a':
			return 'A';
		case 'c':
			return 'C';
		case 'g':
			return 'G';
		case 't':
			return 'T';
		case 'n':
			return 'N';
		default:
			return base;
		}
	}

	/**
	 * Returns the index of the most common base in the basecounts array. To be
	 * used with pileup.getBaseCounts.
	 */
	static public int mostFrequentBaseIndex(int[] baseCounts) {
		int mostFrequentBaseIndex = 0;
		for (int baseIndex = 1; baseIndex < 4; baseIndex++) {
			if (baseCounts[baseIndex] > baseCounts[mostFrequentBaseIndex]) {
				mostFrequentBaseIndex = baseIndex;
			}
		}
		return mostFrequentBaseIndex;
	}

	static public int mostFrequentBaseIndexNotRef(int[] baseCounts,
			int refBaseIndex) {
		int tmp = baseCounts[refBaseIndex];
		baseCounts[refBaseIndex] = -1;
		int result = mostFrequentBaseIndex(baseCounts);
		baseCounts[refBaseIndex] = tmp;
		return result;
	}

	static public int mostFrequentBaseIndexNotRef(int[] baseCounts,
			byte refSimpleBase) {
		return mostFrequentBaseIndexNotRef(baseCounts,
				simpleBaseToBaseIndex(refSimpleBase));
	}

	/**
	 * Returns the most common base in the basecounts array. To be used with
	 * pileup.getBaseCounts.
	 */
	static public byte mostFrequentSimpleBase(int[] baseCounts) {
		return baseIndexToSimpleBase(mostFrequentBaseIndex(baseCounts));
	}

	/**
	 * For the most frequent base in the sequence, return the percentage of the
	 * read it constitutes.
	 */
	static public double mostFrequentBaseFraction(byte[] sequence) {
		int[] baseCounts = new int[4];

		for (byte base : sequence) {
			int baseIndex = simpleBaseToBaseIndex(base);

			if (baseIndex >= 0) {
				baseCounts[baseIndex]++;
			}
		}

		int mostFrequentBaseIndex = mostFrequentBaseIndex(baseCounts);

		return ((double) baseCounts[mostFrequentBaseIndex])
				/ ((double) sequence.length);
	}

	/**
	 * Computes the smallest period >= minPeriod for the specified string. The
	 * period is defined as such p, that for all i = 0... seq.length-1, seq[ i %
	 * p ] = seq[i] (or equivalently seq[i] = seq[i+p] for
	 * i=0...seq.length-1-p). The sequence does <i>not</i> have to contain whole
	 * number of periods. For instance, "ACACACAC" has a period of 2 (it has a
	 * period of 4 as well), and so does "ACACA"; similarly, smallest periods of
	 * "CTCCTC", "CTCCT", and "CTCC" are all equal to 3. The "trivial" period is
	 * the length of the string itself, and it will always be returned if no
	 * smaller period can be found in the specified period range or if specified
	 * minPeriod is greater than the sequence length.
	 *
	 * @param seq
	 * @return
	 */
	public static int sequencePeriod(byte[] seq, int minPeriod) {
		int period = (minPeriod > seq.length ? seq.length : minPeriod);
		// we assume that bases [0,period-1] repeat themselves and check this
		// assumption
		// until we find correct period

		for (int pos = period; pos < seq.length; pos++) {

			int offset = pos % period; // we are currenlty 'offset' bases into
										// the putative repeat of period
										// 'period'
			// if our current hypothesis holds, base[pos] must be the same as
			// base[offset]

			if (Character.toUpperCase(seq[pos]) != Character
					.toUpperCase(seq[offset])) {

				// period we have been trying so far does not work.
				// two possibilities:
				// A) offset = 0, i.e. current position pos must be start of the
				// next repeat, but it is not;
				// in this case only bases from start up to the current one,
				// inclusive, may form a repeat, if at all;
				// so period is at least pos+1 (remember, pos is 0-based), then
				// on the next loop re-entrance
				// pos will be autoincremented and we will be checking next base
				// B) offset != 0, i.e. the current base breaks the repeat, but
				// maybe it starts a new one?
				// hence we should first check if it matches the first base of
				// the sequence, and to do that
				// we set period to pos (thus trying the hypothesis that bases
				// from start up to the current one,
				// non-inclusive are repeated hereafter), and decrement pos
				// (this will re-test current base against the first base
				// on the next loop re-entrance after pos is autoincremented)
				if (offset == 0)
					period = pos + 1;
				else
					period = pos--;

			}
		}
		return period;
	}
}
