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
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect.ErrorWarningType;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;
import org.bgi.flexlab.gaea.tools.annotator.interval.Exon;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant;

/**
 * Calculate codon changes produced by a SNP
 * @author pcingola
 */
public class CodonChangeSnp extends CodonChange {

	public CodonChangeSnp(Variant variant, Transcript transcript, VariantEffects variantEffects) {
		super(variant, transcript, variantEffects);
		returnNow = true; // A SNP can only affect one exon
	}

	/**
	 * Analyze SNPs in this transcript.
	 * Add changeEffect to 'changeEffect'
	 */
	@Override
	protected boolean codonChange(Exon exon) {
		// Get old and new codons
		codonsRef = codonsRef();
		codonsAlt = codonsAlt();

		// Use a generic low priority variant, this allows 'setCodons' to override it
		effect(exon, EffectType.CODON_CHANGE, true);

		if (codonsRef.isEmpty()) variantEffects.addErrorWarning(variant, ErrorWarningType.ERROR_MISSING_CDS_SEQUENCE);

		return true;
	}

	/**
	 * Get new (modified) codons
	 */
	@Override
	protected String codonsAlt() {
		// Was there a problem getting 'codonsOld'? => We cannot do anything
		if (codonsRef.isEmpty()) return "";

		char codonChars[] = codonsRef.toLowerCase().toCharArray();
		char snpBase = variant.netChange(transcript.isStrandMinus()).charAt(0);
		if (codonStartIndex < codonChars.length) codonChars[codonStartIndex] = Character.toUpperCase(snpBase);

		String codonsNew = new String(codonChars);
		return codonsNew;
	}

	/**
	 * Get original codons in CDS
	 */
	@Override
	protected String codonsRef() {
		int numCodons = 1;

		// Get CDS
		String cdsStr = transcript.cds();
		int cdsLen = cdsStr.length();

		// Calculate minBase (first codon base in the CDS)
		int minBase = codonStartNum * CodonChange.CODON_SIZE;
		if (minBase < 0) minBase = 0;

		// Calculate maxBase (last codon base in the CDS)
		int maxBase = codonStartNum * CodonChange.CODON_SIZE + numCodons * CodonChange.CODON_SIZE;
		if (maxBase > cdsLen) maxBase = cdsLen;

		// Sanity checks
		if (cdsStr.isEmpty() // Empty CDS => Cannot get codon (e.g. one or more exons are missing their sequences
				|| (cdsLen <= minBase) // Codon past CDS sequence => Cannot get codon
		) return "";

		// Create codon sequence
		char codonChars[] = cdsStr.substring(minBase, maxBase).toLowerCase().toCharArray();

		// Capitatlize changed base
		if (codonStartIndex < codonChars.length) codonChars[codonStartIndex] = Character.toUpperCase(codonChars[codonStartIndex]);
		String codon = new String(codonChars);

		return codon;
	}
}
