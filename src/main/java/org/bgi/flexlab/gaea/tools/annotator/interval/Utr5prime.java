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
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant.VariantType;
import org.bgi.flexlab.gaea.tools.annotator.util.GprSeq;

import java.util.Collections;
import java.util.List;

/**
 * Interval for a UTR (5 prime UTR and 3 prime UTR
 *
 * @author pcingola
 *
 */
public class Utr5prime extends Utr {

	private static final long serialVersionUID = 3710420226746056364L;
	List<Utr5prime> utrs;

	public Utr5prime() {
		super();
		type = EffectType.UTR_5_PRIME;
	}

	public Utr5prime(Exon parent, int start, int end, boolean strandMinus, String id) {
		super(parent, start, end, strandMinus, id);
		type = EffectType.UTR_5_PRIME;
	}

	synchronized List<Utr5prime> get5primeUtrs() {
		if (utrs == null) {
			Transcript tr = (Transcript) findParent(Transcript.class);

			// Get UTRs and sort them
			utrs = tr.get5primeUtrs();
			if (isStrandPlus()) Collections.sort(utrs, new IntervalComparatorByStart()); // Sort by start position
			else Collections.sort(utrs, new IntervalComparatorByEnd(true)); // Sort by end position (reversed)
		}

		return utrs;
	}

	public String getSequence() {
		// Create UTR sequence
		StringBuffer sb = new StringBuffer();
		for (Utr5prime utr : get5primeUtrs()) {
			Exon ex = (Exon) utr.getParent();
			String utrSeq = ex.getSequence();
			if (utr.size() < utrSeq.length()) utrSeq = utrSeq.substring(0, utr.size()); // UTR5' may stop before end of exon
			sb.append(utrSeq);
		}

		return sb.toString();
	}

	@Override
	public boolean isUtr3prime() {
		return false;
	}

	@Override
	public boolean isUtr5prime() {
		return true;
	}

	/**
	 * Is a new start codon produced?
	 * @param chars
	 * @param pos
	 * @return New start codon (or empty string if there is no new start codon)
	 */
	String startGained(char[] chars, int pos) {
		CodonTable ctable = CodonTables.getInstance().getTable(getGenome(), getChromosomeName());

		// Analyze all frames
		for (int i = Math.max(0, pos - 2); (i <= pos) && ((i + 2) < chars.length); i++) {
			String codon = "" + chars[i] + chars[i + 1] + chars[i + 2];
			if (ctable.isStart(codon)) return codon.toUpperCase(); // This frame has a start codon?
		}
		return "";
	}

	/**
	 * Did we gain a start codon in this 5'UTR interval?
	 * @param seqChange
	 * @return A new start codon (if gained)
	 */
	String startGained(Variant seqChange, Transcript tr) {
		if (!seqChange.isSnp()) return ""; // Only SNPs supported.

		// Calculate SNP position relative to UTRs
		int pos = seqChange.distanceBases(get5primeUtrs(), isStrandMinus());

		// Change base at SNP position
		String sequence = getSequence();
		if(sequence.isEmpty())
			return "";  //TODO
		char[] chars = sequence.toCharArray();
		char snpBase = seqChange.netChange(this).charAt(0);
		if (isStrandMinus()) snpBase = GprSeq.wc(snpBase);
		if(pos >= chars.length) {
			System.err.println("Utr5prime:" + sequence + ";" + pos);
			return "";
		}
		chars[pos] = snpBase;

		// Do we gain a new start codon?
		return startGained(chars, pos);
	}

	/**
	 * Calculate distance from the end of 5'UTRs
	 */
	@Override
	int utrDistance(Variant variant, Transcript tr) {
		int cdsStart = tr.getCdsStart();
		if (cdsStart < 0) return -1;

		if (isStrandPlus()) return cdsStart - variant.getEnd();
		return variant.getStart() - cdsStart;
	}

	@Override
	public boolean variantEffect(Variant variant, VariantEffects variantEffects) {
		// Has the whole UTR been deleted?
		if (variant.includes(this) && (variant.getVariantType() == VariantType.DEL)) {
			variantEffects.add(variant, this, EffectType.UTR_5_DELETED, ""); // A UTR was removed entirely
			return true;
		}

		// Add distance
		Transcript tr = (Transcript) findParent(Transcript.class);
		int distance = utrDistance(variant, tr);
		VariantEffect variantEffect = new VariantEffect(variant);
		variantEffect.set(this, type, type.effectImpact(), distance >= 0 ? distance + " bases from TSS" : "");
		variantEffect.setDistance(distance);
		variantEffects.add(variantEffect);

		// Start gained?
		String gained = startGained(variant, tr);
		if (!gained.isEmpty()) variantEffects.add(variant, this, EffectType.START_GAINED, gained);

		return true;
	}

}
