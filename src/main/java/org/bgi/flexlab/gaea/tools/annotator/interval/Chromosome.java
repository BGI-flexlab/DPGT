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

import org.bgi.flexlab.gaea.tools.annotator.codons.CodonTable;
import org.bgi.flexlab.gaea.tools.annotator.codons.CodonTables;
import org.bgi.flexlab.gaea.tools.annotator.effect.EffectType;
import org.bgi.flexlab.gaea.tools.annotator.util.Gpr;

/**
 * Interval for the whole chromosome
 * If a SNP has no 'ChromosomeInterval' => it is outside the chromosome => Invalid
 *
 * @author pcingola
 *
 */
public class Chromosome extends Marker {

	private static final long serialVersionUID = 1636197649250882952L;

	double chromosomeNum;
//	DnaSequence sequence = null;

	/**
	 * Compare chromosome names
	 */
	public static int compare(String chr1, String chr2) {
		// Try to compare numbers
		int chr1Num = number(chr1);
		int chr2Num = number(chr2);
		if (chr1Num > 0 && chr2Num > 0) return chr1Num - chr2Num;
		if (chr1Num > 0 && chr2Num <= 0) return -1;
		if (chr1Num <= 0 && chr2Num > 0) return 1;

		// Numbers did not work, compare strings
		return simpleName(chr1).compareTo(simpleName(chr2));
	}

	/**
	 * Convert to chromosome number (return '0' if it cannot be converted)
	 */
	public static int number(String chrName) {
		return Gpr.parseIntSafe(ChromosomeSimpleName.get(chrName));
	}

	/**
	 * Simplify chromosome name
	 */
	public static String simpleName(String chrName) {
		return ChromosomeSimpleName.get(chrName);
	}

	public Chromosome() {
		super();
	}

	public Chromosome(Genome parent, int start, int end, String id) {
		super(null, start, end, false, id); // Parent = null to avoid sanity check (it will always fail for chromosomes)
		this.parent = parent;
		type = EffectType.CHROMOSOME;
		setChromosomeName(id);
	}

	@Override
	public Chromosome cloneShallow() {
		Chromosome clone = (Chromosome) super.cloneShallow();
		clone.chromosomeNum = chromosomeNum;
		return clone;
	}

	/**
	 * Compare only chromosome's name
	 */
	public int compareChromoName(Interval interval) {
		Chromosome i2 = (Chromosome) interval;

		// Both of them are non-numeric: Compare as string
		if ((chromosomeNum == 0) && (i2.chromosomeNum == 0)) return id.compareTo(i2.id);

		// One of them is a number? => the number goes first
		if ((chromosomeNum != 0) && (i2.chromosomeNum == 0)) return -1;
		if ((chromosomeNum == 0) && (i2.chromosomeNum != 0)) return 1;

		// Use numeric comparison
		if (chromosomeNum - i2.chromosomeNum < 0) return -1;
		if (chromosomeNum - i2.chromosomeNum > 0) return 1;
		return 0;
	}

	public CodonTable getCodonTable() {
		return CodonTables.getInstance().getTable(getGenome(), getId());
	}

//	public DnaSequence getDnaSequence() {
//		return sequence;
//	}
//
//	public String getSequence() {
//		return sequence.toString();
//	}

	/**
	 * Is this a mitochondrial chromosome?
	 * Note: This is a wild guess just by looking at the name
	 */
	public boolean isMt() {
		String iduc = id.toUpperCase();
		return iduc.equals("M") //
				|| iduc.startsWith("MT") //
				|| (iduc.indexOf("MITO") >= 0) //
				;
	}

	/**
	 * Set chromosome name
	 * Note: Removes prefixes (such as 'chr') and parse numeric version.
	 * E.g. 'chr2' becomes '2'. Also numeric '2' is assigned to 'chromosomeNum' to facilitate order by number (so that '2' is ordered before '21')
	 * @param chromo
	 */
	private void setChromosomeName(String chromo) {
		id = simpleName(chromo);
		chromosomeNum = Gpr.parseIntSafe(id); // Try to parse a numeric string
	}

	public void setLength(int len) {
		end = len - 1; // Remember that intervals are zero-based
	}

//	/**
//	 * Set sequence for this chromosome
//	 * @param sequenceStr
//	 */
//	public void setSequence(String sequenceStr) {
//		sequence = new DnaSequence(sequenceStr, true);
//		setLength(sequenceStr.length()); // Update chromosome length
//	}

}
