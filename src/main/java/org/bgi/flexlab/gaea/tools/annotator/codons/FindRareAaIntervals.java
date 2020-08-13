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
package org.bgi.flexlab.gaea.tools.annotator.codons;

import org.bgi.flexlab.gaea.tools.annotator.interval.Gene;
import org.bgi.flexlab.gaea.tools.annotator.interval.Genome;
import org.bgi.flexlab.gaea.tools.annotator.interval.RareAminoAcid;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.util.GprSeq;

import java.util.HashMap;

/**
 * Find intervals where rare amino acids occur
 * 
 * @author pablocingolani
 */
public class FindRareAaIntervals {

	public static final double RARE_THRESHOLD = 1e-5;

	boolean verbose = false;
	double rareThreshold = RARE_THRESHOLD;
	Genome genome;
	CodonTable codonTable;
	int count[];
	int countTotal;
	boolean isInTable[];
	HashMap<String, Transcript> trById;
	HashMap<String, RareAminoAcid> rareAaByPos;

	public FindRareAaIntervals(Genome genome) {
		this.genome = genome;
		codonTable = genome.codonTable();
		count = new int['Z']; // It's all upper case
		isInTable = new boolean['Z']; // It's all upper case
		countTotal = 0;
		rareAaByPos = new HashMap<String, RareAminoAcid>();
	}

	/**
	 * Add a rare amino acid (if not already added)
	 * @param tr
	 * @param start
	 * @param end
	 */
	void addRareAa(Transcript tr, int start, int end) {
		int s = Math.min(start, end);
		int e = Math.max(start, end);

		String key = tr.getChromosomeName() + ":" + s + "-" + e;
		if( rareAaByPos.containsKey(key) ) return; // Nothing to add

		RareAminoAcid raa = new RareAminoAcid(tr, s, e, "");
		rareAaByPos.put(key, raa);
	}

	/**
	 * Find interval AA index in a transcript
	 * @param trId
	 * @param aaIdx
	 */
	void findRareAaSites(String trId, int aaIdx) {
		// Create index?
		if( trById == null ) {
			trById = new HashMap<String, Transcript>();
			for( Gene gene : genome.getGenome().getGenes() )
				for( Transcript tr : gene )
					trById.put(tr.getId(), tr);
		}

		// Find transcript
		Transcript tr = trById.get(trId);
		if( tr == null ) {
			if( verbose ) System.err.println("WARNING: Cannot find transcript '" + trId + "'");
			return;
		}

		// Create markers. There might be more than one interval since an intron can be in the middle of the codon.
		int cds2pos[] = tr.baseNumberCds2Pos();
		int pos = 0, posPrev = 0, start = -1;
		int step = tr.isStrandPlus() ? 1 : -1;
		for( int cds = aaIdx * 3; cds < (aaIdx + 1) * 3; cds++ ) {
			pos = cds2pos[cds];
			if( start < 0 ) start = pos;
			else if( pos != posPrev + step ) {
				// Non-contiguous: Create a new marker and add it to the list
				addRareAa(tr, start, pos);
				start = -1;
			}

			posPrev = pos;
		}

		if( start >= 0 ) addRareAa(tr, start, pos);
	}

	/**
	 * A string containing all rare amino acids
	 * @return
	 */
	String findRareNames() {
		StringBuilder rare = new StringBuilder();
		for( int i = 0; i < count.length; i++ ) {
			double p = ((double) count[i]) / countTotal;
			if( (count[i] > 0) && ((p < rareThreshold) || !isInTable[i]) ) rare.append((char) i);
		}
		return rare.toString();
	}

	public double getRareThreshold() {
		return rareThreshold;
	}

	/**
	 * Calculate all AA in a codon table
	 */
	void isInTable() {
		for( char c1 : GprSeq.BASES )
			for( char c2 : GprSeq.BASES )
				for( char c3 : GprSeq.BASES ) {
					String codon = "" + c1 + c2 + c3;
					String aa = codonTable.aa(codon);
					isInTable[aa.toUpperCase().charAt(0)] = true;
				}
	}

	public void setRareThreshold(double rareThreshold) {
		this.rareThreshold = rareThreshold;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for( int i = 0; i < count.length; i++ ) {
			double p = ((double) count[i]) / countTotal;
			if( count[i] > 0 ) sb.append(String.format("\t%s\t%d\t%.2e\t%b\n", (char) i, count[i], p, isInTable[i]));
		}

		return sb.toString();
	}
}
