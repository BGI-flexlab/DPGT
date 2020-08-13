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
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect.EffectImpact;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;
import org.bgi.flexlab.gaea.tools.annotator.interval.Exon;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant;

/**
 * Calculate codon changes produced by a deletion
 * @author pcingola
 */
public class CodonChangeDel extends CodonChangeStructural {

	public CodonChangeDel(Variant variant, Transcript transcript, VariantEffects variantEffects) {
		super(variant, transcript, variantEffects);
		returnNow = false;
		requireNetCdsChange = true;
	}

	/**
	 * Analyze deletions in this transcript.
	 */
	@Override
	protected boolean codonChange(Exon exon) {
		// Is there any net effect?
		if (netCdsChange.isEmpty()) return false;

		EffectType effType = null;

		if (variant.includes(exon)) {
			/**
			 * An exon has been entirely removed
			 */
			codonsRef = "";
			codonsAlt = "";
			codonStartNum = codonStartIndex = -1;
			effType = EffectType.EXON_DELETED;
		} else if (netCdsChange.length() % CodonChange.CODON_SIZE != 0) {
			/**
			 * Length not multiple of CODON_SIZE => FRAME_SHIFT
			 * 	E.g. :
			 * 		Original:			AAA CCC GGG AAA CCC GGG AAA CCC GGG
			 * 		Delete 'AA' pos 0:	ACC CGG GAA ACC CGG GAA ACC CGG G
			 * 		Delete 'AA' pos 1:	ACC CGG GAA ACC CGG GAA ACC CGG G
			 * 		Delete 'AC' pos 2:	AAC CGG GAA ACC CGG GAA ACC CGG G
			 */
			codonsRef = codonsRef();
			codonsAlt = "";
			effType = EffectType.FRAME_SHIFT;
		} else if (codonStartIndex == 0) {
			/**
			 * Length multiple of CODON_SIZE and insertion happens at codon boundary => CODON_INSERTION
			 * 	E.g. :
			 * 		Original:			AAA CCC GGG AAA CCC GGG AAA CCC GGG
			 * 		Delete 'AAA' pos 0:	CCC GGG AAA CCC GGG AAA CCC GGG
			 */
			codonsRef = codonsRef();
			codonsAlt = "";
			effType = EffectType.CODON_DELETION;
		} else {
			/**
			 * Length multiple of CODON_SIZE and insertion does not happen at codon boundary => CODON_CHANGE_PLUS_CODON_DELETION
			 * 	E.g. :
			 * 		Original:			AAA CCC GGG AAA CCC GGG AAA CCC GGG
			 * 		Delete 'AAC' pos 1:	ACC GGG AAA CCC GGG AAA CCC GGG
			 * 		Delete 'ACC' pos 2:	AAC GGG AAA CCC GGG AAA CCC GGG
			 */
			codonsRef = codonsRef();
			codonsAlt = codonsAlt();

			if (codonsAlt.isEmpty() || codonsRef.startsWith(codonsAlt)) {
				/**
				 * Note: It might happen that the last codon of the exon was deleted.
				 *       In this case there is no 'CODON_CHANGE'
				 * E.g.
				 * 		Original:				AAA CCC GGG AAA CCC GGG AAA CCC GGG
				 * 		Delete 'GGG' pos 24:	ACC CCC GGG AAA CCC GGG AAA CCC
				 *
				 * Note2: It may also be the case that the deleted bases are equal to the following ones.
				 *  E.g.
				 *  	Original:			ACG TCG TCC GGG AAA CCC GGG AAA CCC GGG
				 *  	Delete 'CGT' pos 1:	ACG TCC GGG AAA CCC GGG AAA CCC GGG
				 */
				effType = EffectType.CODON_DELETION;
			} else {
				effType = EffectType.CODON_CHANGE_PLUS_CODON_DELETION;
			}
		}

		effect(exon, effType, false);

		return true;
	}

	/**
	 * Get alternative codons
	 */
	@Override
	protected String codonsAlt() {
		if (netCdsChange.isEmpty()) return "";

		int after = netCdsChange.length() + codonStartIndex;

		String prefix = codonsRef.length() >= codonStartIndex ? codonsRef.substring(0, codonStartIndex) : codonsRef;
		String suffix = codonsRef.length() > after ? codonsRef.substring(after) : "";

		String codonsAlt = prefix + suffix;
		return codonsAlt;
	}

	/**
	 * Get original codons in CDS
	 */
	@Override
	protected String codonsRef() {
		if (netCdsChange.isEmpty()) return "";

		int min = variant.getStart();
		int max = variant.getEnd();
		int cdsBaseMin = cdsBaseNumber(min);
		int cdsBaseMax = cdsBaseNumber(max);

		// Swap?
		if (transcript.isStrandMinus()) {
			int swap = cdsBaseMin;
			cdsBaseMin = cdsBaseMax;
			cdsBaseMax = swap;
		}

		if (cdsBaseMax < cdsBaseMin) throw new RuntimeException("This should never happen!\n\tcdsBaseMin: " + cdsBaseMin + "\n\tcdsBaseMax: " + cdsBaseMax + "\n\tmin: " + min + "\n\tmax: " + max + "\n\tSeqChange: " + variant + "\n\ttranscript: " + transcript + "\n\tCDS.len: " + transcript.cds().length());

		int maxCodon = cdsBaseMax / CodonChange.CODON_SIZE;
		int minCodon = cdsBaseMin / CodonChange.CODON_SIZE;
		int oldCodonCdsStart = (CodonChange.CODON_SIZE * minCodon);
		int oldCodonCdsEnd = (CodonChange.CODON_SIZE * (maxCodon + 1)) - 1;

		String codons = "";
		if (oldCodonCdsEnd >= transcript.cds().length()) codons = transcript.cds().substring(oldCodonCdsStart);
		else codons = transcript.cds().substring(oldCodonCdsStart, oldCodonCdsEnd + 1);

		return codons;
	}

	@Override
	protected void effectTranscript() {
		effectNoCodon(transcript, EffectType.TRANSCRIPT_DELETED);
	}

	/**
	 * Whole exon/s deleted?
	 */
	void exonLoss() {
		for (Exon ex : transcript)
			if (variant.includes(ex)) effectNoCodon(ex, EffectType.EXON_DELETED);
	}

	/**
	 * Deletion analysis using full transcript information. This is
	 * done only when the variant affects more than one exons.
	 */
	@Override
	protected void exons() {
		if (exonFull == 0 && exonPartial == 1) {
			// Variant partially affects only one exon?
			// => Use the standard (by exon) method
			codonChangeSuper();
			return;
		} else if (exonFull > 0) {
			// Full exons deleted
			exonLoss();

			// Only whole exons deleted? We are done
			if (exonPartial == 0) return;
		}

		//---
		// A combination of partial and full exons affected
		//---
		codonsRefAlt();
		EffectType effType = null;
		int lenDiff = cdsAlt.length() - cdsRef.length();
		if (lenDiff % CodonChange.CODON_SIZE != 0) {
			effType = EffectType.FRAME_SHIFT;
		} else if (codonStartIndex == 0) {
			effType = EffectType.CODON_DELETION;
		} else {
			if (codonsAlt.isEmpty() || codonsRef.startsWith(codonsAlt)) {
				effType = EffectType.CODON_DELETION;
			} else {
				effType = EffectType.CODON_CHANGE_PLUS_CODON_DELETION;
			}
		}

		// Assign to first exon
		for (Exon ex : transcript)
			if (variant.includes(ex) || variant.intersects(ex)) {
				exon = ex;
				break;
			}

		// Add variant effect
		effect(exon, effType, false);
	}

	@Override
	protected void exonsCoding() {
		// TODO Auto-generated method stub

	}

	@Override
	protected void exonsNoncoding() {
		if (exonFull > 0) effectNoCodon(transcript, EffectType.EXON_DELETED, EffectImpact.MODIFIER);
		if (exonPartial > 0) effectNoCodon(transcript, EffectType.EXON_DELETED_PARTIAL, EffectImpact.MODIFIER);
	}

	@Override
	protected void intron() {
		effectNoCodon(transcript, EffectType.INTRON);
	}

}
