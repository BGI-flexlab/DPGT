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
package org.bgi.flexlab.gaea.tools.annotator.interval.codonchange;

import org.bgi.flexlab.gaea.tools.annotator.config.Config;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;
import org.bgi.flexlab.gaea.tools.annotator.interval.Exon;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant;
import org.bgi.flexlab.gaea.tools.annotator.util.Gpr;

/**
 * Calculate codon changes produced by a duplication
 *
 * @author pcingola
 */
public abstract class CodonChangeStructural extends CodonChange {

	public static boolean debug = false;

	protected boolean coding;
	protected int exonFull, exonPartial;
	protected String cdsAlt;
	protected String cdsRef;

	public CodonChangeStructural(Variant variant, Transcript transcript, VariantEffects variantEffects) {
		super(variant, transcript, variantEffects);
		coding = transcript.isProteinCoding() || Config.get().isTreatAllAsProteinCoding();
		countAffectedExons();
	}

	/**
	 * Differences between two CDSs after removing equal codons from
	 * the beginning and from the end of both strings
	 */
	protected void cdsDiff() {
		int min = Math.min(cdsRef.length(), cdsAlt.length()) / 3;

		// Removing form the beginning
		codonStartNum = 0;
		codonStartIndex = 0;
		for (int i = 0; i <= min; i++) {
			codonStartNum = i;
			if (debug) Gpr.debug("cdsDiff Start\tcodonEquals(" + i + " , " + i + "): " + codonEquals(cdsRef, cdsAlt, i, i) //
					+ "\n\tcodonsRef [" + codonStartNum + "]: " + codons(cdsRef, codonStartNum, -1) //
					+ "\n\tcodonsAlt [" + codonStartNum + "]: " + codons(cdsAlt, codonStartNum, -1) //
			);

			// Codons differ?
			if (!codonEquals(cdsRef, cdsAlt, i, i)) {
				// Find index difference within codon
				codonStartIndex = codonDiffIndex(cdsRef, cdsAlt, i, i);
				break;
			}
		}

		// Removing trailing codons
		int codonNumEndRef = cdsRef.length() / 3; //+ (cdsRef.length() % 3 == 0 ? 0 : 1);
		int codonNumEndAlt = cdsAlt.length() / 3; //+ (cdsAlt.length() % 3 == 0 ? 0 : 1);

		for (; codonNumEndRef >= codonStartNum && codonNumEndAlt >= codonStartNum; codonNumEndRef--, codonNumEndAlt--) {
			if (debug) Gpr.debug("cdsDiff End\tcodonEquals(" + codonNumEndRef + " , " + codonNumEndAlt + "): " + codonEquals(cdsRef, cdsAlt, codonNumEndRef, codonNumEndAlt) //
					+ "\n\tcodonsRef [" + codonStartNum + " , " + codonNumEndRef + "]: " + codons(cdsRef, codonStartNum, codonNumEndRef) //
					+ "\n\tcodonsAlt [" + codonStartNum + " , " + codonNumEndAlt + "]: " + codons(cdsAlt, codonStartNum, codonNumEndAlt) //
			);
			if (!codonEquals(cdsRef, cdsAlt, codonNumEndRef, codonNumEndAlt)) break;
		}

		// Codons
		codonsRef = codons(cdsRef, codonStartNum, codonNumEndRef);
		codonsAlt = codons(cdsAlt, codonStartNum, codonNumEndAlt);

		// No codon difference found?
		if (codonsRef.isEmpty() && codonsAlt.isEmpty()) codonStartNum = codonStartIndex = -1;
	}

	@Override
	public void codonChange() {
		if (variant.includes(transcript)) {
			// Whole transcript affected?
			effectTranscript();
		} else {
			// Does the variant affect any exons?
			if ((exonFull > 0 || exonPartial > 0)) exons();
			else intron();
		}
	}

	protected void codonChangeSuper() {
		super.codonChange();
	}

	private int codonDiffIndex(String cdsRef, String cdsAlt, int codonNumRef, int codonNumAlt) {
		for (int h = 0, i = 3 * codonNumRef, j = 3 * codonNumAlt; h < 3; i++, j++, h++) {
			// Premature end of sequence? (i.e. sequence ends before codon end)
			if ((i >= cdsRef.length()) || (j >= cdsAlt.length())) return h;

			// Different base? Return index
			if (cdsRef.charAt(i) != cdsAlt.charAt(j)) return h;
		}

		return -1;
	}

	/**
	 * Compare codons from cdsRef[codonNumRef] and cdsAlt[codonNumAlt]
	 */
	private boolean codonEquals(String cdsRef, String cdsAlt, int codonNumRef, int codonNumAlt) {
		for (int h = 0, i = 3 * codonNumRef, j = 3 * codonNumAlt; h < 3; i++, j++, h++) {
			if ((i >= cdsRef.length()) || (j >= cdsAlt.length())) // Premature end of sequence? (i.e. sequence ends before codon end)
				return (i >= cdsRef.length()) && (j >= cdsAlt.length()); // We consider them equal only if both sequences reached the end at the same time

			// Same base?
			if (cdsRef.charAt(i) != cdsAlt.charAt(j)) return false;
		}

		return true;
	}

	/**
	 * Get codons from CDS
	 */
	private String codons(String cds, int codonNumStart, int codonNumEnd) {
		if (codonNumEnd >= 0 && codonNumEnd < codonNumStart) return "";
		int endBase = codonNumEnd < 0 ? cds.length() : 3 * (codonNumEnd + 1);
		endBase = Math.min(cds.length(), endBase);
		int startBase = Math.max(0, 3 * codonNumStart);
		return cds.substring(startBase, endBase);
	}

	/**
	 * Calculate codons by applying the variant and calculating the differences in CDS sequences
	 * This is a slow method, makes sense only for complex variants
	 */
	protected void codonsRefAlt() {
		Transcript trNew = transcript.apply(variant);
		if (debug) Gpr.debug("Transcript after apply: " + trNew);

		cdsAlt = trNew.cds();
		cdsRef = transcript.cds();

		// Calculate differences: CDS
		cdsDiff();
	}

	/**
	 * How many full / partial exons does the variant affect?
	 */
	protected void countAffectedExons() {
		exonFull = exonPartial = 0;

		for (Exon ex : transcript) {
			if (variant.includes(ex)) exonFull++;
			else if (variant.intersects(ex)) exonPartial++;
		}
	}

	protected abstract void effectTranscript();

	/**
	 * Variant affect one or more exons
	 */
	protected abstract void exons();

	/**
	 * Variant affect one or more coding exons
	 */
	protected abstract void exonsCoding();

	/**
	 * Variant affect one or more non-coding exons
	 */
	protected abstract void exonsNoncoding();

	/**
	 * Variant affect one intron
	 */
	protected abstract void intron();

}
