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
 * Copyright (C)  2016  Pablo Cingolani(pcingola@users.sourceforge.net)
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.annotator.util;

import java.util.Iterator;

/**
 * Find all bases combinations from a string containing IUB codes
 *
 * @author pcingola
 */
public class IubString implements Iterable<String>, Iterator<String> {

	// Maximum number of IUB bases allowed
	public static final int MAX_IUB_BASES = 10;

	int idx2BaseNum[];
	char bases[];
	char iubCodesByIndex[][];

	CombinatorialIterator combIt;

	/**
	 * How many IUB bases are in this string?
	 */
	public static int countIUB(String str) {
		char bases[] = str.toCharArray();

		int count = 0;
		for (int i = 0; i < bases.length; i++)
			if (isUIB(bases[i])) return count++;

		return count;
	}

	/**
	 * Does the string have ANY IUB base?
	 */
	public static boolean hasIUB(String str) {
		char bases[] = str.toCharArray();

		for (int i = 0; i < bases.length; i++)
			if (isUIB(bases[i])) return true;

		return false;
	}

	/**
	 * Does the string have at most 'MAX_IUB_BASES' IUB bases?
	 */
	public static boolean hasIUBMax(String str) {
		if (str.length() < MAX_IUB_BASES) return hasIUB(str);
		int count = countIUB(str);
		return (count > 0) && (count <= MAX_IUB_BASES);
	}

	public static boolean isUIB(char base) {
		switch (base) {
		case 'N': // aNy base
		case 'B': // B: not A
		case 'D': // D: not C
		case 'H': // H: not G
		case 'V': // V: not T
		case 'M':
		case 'R':
		case 'W': // Weak
		case 'S': // Strong
		case 'Y':
		case 'K':
			return true;

		default:
			return false;
		}
	}

	/**
	 * Convert a single IUB code to the corresponding bases
	 *  IUB codes: M=A/C, R=A/G, W=A/T, S=C/G, Y=C/T, K=G/T and N=A/C/G/T
	 */
	public static char[] iub2bases(char alt) {
		char[] alts;

		switch (alt) {
		case 'N': // aNy base
			alts = new char[4];
			alts[0] = 'A';
			alts[1] = 'C';
			alts[2] = 'G';
			alts[3] = 'T';
			break;

		case 'B': // B: not A
			alts = new char[3];
			alts[0] = 'C';
			alts[1] = 'G';
			alts[2] = 'T';
			break;

		case 'D': // D: not C
			alts = new char[3];
			alts[0] = 'A';
			alts[1] = 'G';
			alts[2] = 'T';
			break;

		case 'H': // H: not G
			alts = new char[3];
			alts[0] = 'A';
			alts[1] = 'C';
			alts[2] = 'T';
			break;

		case 'V': // V: not T
			alts = new char[3];
			alts[0] = 'A';
			alts[1] = 'C';
			alts[2] = 'G';
			break;

		case 'M':
			alts = new char[2];
			alts[0] = 'A';
			alts[1] = 'C';
			break;

		case 'R':
			alts = new char[2];
			alts[0] = 'A';
			alts[1] = 'G';
			break;

		case 'W': // Weak
			alts = new char[2];
			alts[0] = 'A';
			alts[1] = 'T';
			break;

		case 'S': // Strong
			alts = new char[2];
			alts[0] = 'C';
			alts[1] = 'G';
			break;

		case 'Y':
			alts = new char[2];
			alts[0] = 'C';
			alts[1] = 'T';
			break;

		case 'K':
			alts = new char[2];
			alts[0] = 'G';
			alts[1] = 'T';
			break;

		default:
			throw new RuntimeException("WARNING: Unkown IUB code for SNP '" + alt + "'");
		}

		return alts;
	}

	public IubString(String str) {
		bases = str.toCharArray();

		// Get combinatorial iterator size
		int size = 0;
		for (int i = 0; i < bases.length; i++)
			if (isUIB(bases[i])) size++;

		// Initialize
		combIt = new CombinatorialIterator(size);
		iubCodesByIndex = new char[size][];
		idx2BaseNum = new int[size];

		// Initialize iterator parameters
		int j = 0;
		for (int i = 0; i < bases.length; i++) {
			if (isUIB(bases[i])) {
				// Transform base to codes
				char iubCodes[] = iub2bases(bases[i]);
				combIt.set(j, 0, iubCodes.length - 1);

				// Store codes by counter index
				iubCodesByIndex[j] = iubCodes;

				// Map counter's index to base number
				idx2BaseNum[j] = i;
				j++;
			}
		}
	}

	@Override
	public boolean hasNext() {
		return combIt.hasNext();
	}

	@Override
	public Iterator<String> iterator() {
		return this;
	}

	@Override
	public String next() {
		int next[] = combIt.next();
		if (next == null) return null;

		for (int j = 0; j < idx2BaseNum.length; j++) {
			int i = idx2BaseNum[j];
			char iubs[] = iubCodesByIndex[j];
			bases[i] = iubs[next[j]];
		}

		return new String(bases);
	}

	@Override
	public void remove() {
	}

}
