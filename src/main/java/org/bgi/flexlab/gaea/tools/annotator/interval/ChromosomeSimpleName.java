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
package org.bgi.flexlab.gaea.tools.annotator.interval;

import java.util.HashMap;

/**
 * Convert chromosome names to simple names
 * @author pcingola
 */
public class ChromosomeSimpleName {

	public static final String CHROMO_PREFIX[] = { "chromosome", "chromo", "chr" }; //, "group", "scaffold", "contig", "supercontig", "supercont", "0" }; // Must be lower case (see method)
	private static ChromosomeSimpleName instance = new ChromosomeSimpleName();

	private final HashMap<String, String> map;

	/**
	 * Get a simple name for the chromosome
	 */
	public static String get(String chrName) {
		return instance.simpleNameCache(chrName);
	}

	private ChromosomeSimpleName() {
		map = new HashMap<String, String>();
	}

	/**
	 * Remove 'prefix' from 'chr'
	 */
	String removePrefix(String chr, String prefix) {
		String chrLower = chr.toLowerCase();
		String prefixLower = prefix.toLowerCase();

		// Remove prefix it it matches
		if (chrLower.startsWith(prefixLower)) return chr.substring(prefix.length());
		return chr;
	}

	/**
	 * Simplify chromosome name
	 */
	protected String simpleName(String chr) {
		if (chr == null) return "";
		chr = chr.trim();

		// Remove any prefix string until no change is made
		String chrPrev = "";
		do {
			chrPrev = chr;

			// If chr matches EXACTLY, we don't remove any prefix
			// E.g.:
			//		Let's assume that chr is 'Chromosome'
			//			i) If we remove prefix 'Chromosome', then we get an empty chromosome name, which is illegal.
			//			ii) If we remove "Chromo" (another prefix), the we get "some", which doesn't make any sense.
			//
			// So in case of ANY exact match, it's better not to remove any prefix.
			boolean exactMatch = false;
			for (String prefix : CHROMO_PREFIX)
				exactMatch |= chr.equalsIgnoreCase(prefix);

			// Remove all common prefixes
			if (!exactMatch) {
				for (String prefix : CHROMO_PREFIX) {
					chr = removePrefix(chr, prefix + ":");
					chr = removePrefix(chr, prefix + "_");
					chr = removePrefix(chr, prefix + "-");
					chr = removePrefix(chr, prefix);
				}
			}

			// Remove left zeros (unless the name is '0') and trim spaces
			if (!chr.equals("0")) chr = removePrefix(chr, "0");
			chr = chr.trim();

		} while (!chr.equals(chrPrev));

		return chr;
	}

	/**
	 * Query cache before simplifying name
	 */
	protected String simpleNameCache(String chrName) {
		String chr = map.get(chrName);
		if (chr == null) {
			chr = simpleName(chrName);
			map.put(chrName, chr);
		}
		return chr;
	}

}
