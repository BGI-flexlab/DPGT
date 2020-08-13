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

import org.bgi.flexlab.gaea.tools.annotator.effect.EffectType;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;
import org.bgi.flexlab.gaea.tools.annotator.interval.Exon;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant;
import org.bgi.flexlab.gaea.tools.annotator.util.GprSeq;

/**
 * Calculate codon changes produced by a MNP
 * @author pcingola
 */
public class CodonChangeMnp extends CodonChange {

	int cdsStart;
	int cdsEnd;

	public CodonChangeMnp(Variant variant, Transcript transcript, VariantEffects variantEffects) {
		super(variant, transcript, variantEffects);
		returnNow = false;
		requireNetCdsChange = true;
	}

	int cdsBaseNumber(int pos, boolean usePrevBaseIntron) {
		if (pos < cdsStart) return transcript.isStrandPlus() ? 0 : transcript.cds().length() - 1;
		if (pos > cdsEnd) return transcript.isStrandPlus() ? transcript.cds().length() - 1 : 0;
		return transcript.baseNumberCds(pos, usePrevBaseIntron);
	}

	/**
	 * Calculate a list of codon changes
	 */
	@Override
	public void codonChange() {
		codonOldNew();

		// Create change effect
		effect(transcript, EffectType.CODON_CHANGE, true); // Use a generic low priority variant, this allows 'setCodons' to override it

		return;
	}

	/**
	 * Calculate codons old / codons new
	 */
	protected void codonOldNew() {
		if (!transcript.intersects(variant)) return;

		// CDS coordinates
		cdsStart = transcript.isStrandPlus() ? transcript.getCdsStart() : transcript.getCdsEnd();
		cdsEnd = transcript.isStrandPlus() ? transcript.getCdsEnd() : transcript.getCdsStart();

		// Does it intersect CDS?
		if (cdsStart > variant.getEnd()) return;
		if (cdsEnd < variant.getStart()) return;

		// Base number relative to CDS start
		int scStart, scEnd;
		if (transcript.isStrandPlus()) {
			scStart = cdsBaseNumber(variant.getStart(), false);
			scEnd = cdsBaseNumber(variant.getEnd(), true);
		} else {
			scEnd = cdsBaseNumber(variant.getStart(), true);
			scStart = cdsBaseNumber(variant.getEnd(), false);
		}

		// Update coordinates
		codonStartNum = scStart / CODON_SIZE;
		codonStartIndex = scStart % CODON_SIZE;

		// MNP overlap in coding part
		int scLen = scEnd - scStart;
		if (scLen < 0) return;

		// Round to codon position
		int scStart3 = round3(scStart, false);
		int scEnd3 = round3(scEnd, true);
		if (scEnd3 == scStart3) scEnd3 += 3; // At least one codon

		// Append 'N'
		String padN = "";
		int diff = scEnd3 - (transcript.cds().length() - 1);
		if (diff > 0) {
			scEnd3 = transcript.cds().length() - 1;
			// Pad with 'N'
			switch (diff) {
			case 1:
				padN = "N";
				break;

			case 2:
				padN = "NN";
				break;

			default:
				throw new RuntimeException("Sanity check failed. Number of 'N' pading is :" + diff + ". This should not happen!");
			}
		}

		// Get old codon (reference)
		codonsRef = transcript.cds().substring(scStart3, scEnd3 + 1);

		// Get new codon (change)
		String prepend = codonsRef.substring(0, scStart - scStart3);
		String append = "";
		if (scEnd3 > scEnd) append = codonsRef.substring(codonsRef.length() - (scEnd3 - scEnd));
		codonsAlt = prepend + netCdsChange() + append;

		// Pad codons with 'N' if required
		codonsRef += padN;
		codonsAlt += padN;

		//---
		// Can we simplify codons?
		//---
		if ((codonsRef != null) && (codonsAlt != null)) {
			while ((codonsRef.length() >= 3) && (codonsAlt.length() >= 3)) {
				// First codon
				String cold = codonsRef.substring(0, 3);
				String cnew = codonsAlt.substring(0, 3);

				// Are codons equal? => Simplify
				if (cold.equalsIgnoreCase(cnew)) {
					codonsRef = codonsRef.substring(3);
					codonsAlt = codonsAlt.substring(3);
					codonStartNum++;
				} else break;
			}
		}
	}

	@Override
	protected String codonsAlt() {
		return codonsAlt;
	}

	/**
	 * Calculate old codons
	 */
	@Override
	protected String codonsRef() {
		return codonsRef;
	}

	/**
	 * We may have to calculate 'netCdsChange', which is the effect on the CDS
	 * Note: A deletion or a MNP might affect several exons
	 * @return
	 */
	@Override
	protected String netCdsChange() {
		if (variant.size() > 1) {
			StringBuilder sb = new StringBuilder();
			for (Exon exon : transcript.sortedStrand()) {
				String seq = variant.netChange(exon);
				sb.append(exon.isStrandPlus() ? seq : GprSeq.reverseWc(seq));
			}
			return sb.toString();
		}

		return variant.netChange(transcript.isStrandMinus());
	}

	int round3(int num, boolean end) {
		if (end) {
			if (num % 3 == 2) return num;
			return (num / 3) * 3 + 2;
		}

		if (num % 3 == 0) return num;
		return (num / 3) * 3;
	}

}
