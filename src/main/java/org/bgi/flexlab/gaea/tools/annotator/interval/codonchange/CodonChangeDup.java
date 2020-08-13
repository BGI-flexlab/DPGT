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
import org.bgi.flexlab.gaea.tools.annotator.effect.EffectType;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect.EffectImpact;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant;

/**
 * Calculate codon changes produced by a duplication
 *
 * @author pcingola
 */
public class CodonChangeDup extends CodonChangeStructural {

	public CodonChangeDup(Variant variant, Transcript transcript, VariantEffects variantEffects) {
		super(variant, transcript, variantEffects);
		coding = transcript.isProteinCoding() || Config.get().isTreatAllAsProteinCoding();
	}

	/**
	 * Analyze whether the duplication is past transcript's coding region
	 *
	 * E.g.:   Transcript  = chr1:100-200
	 *         Duplication = chr1:150-999
	 *         The duplicated segment starts at base 1000, which is beyond's
	 *         transcript's end, so it probably has no effect on the amino
	 *         acid sequence
	 *
	 *  Rationale:
	 *  If we have two genes:
	 *
	 *     | gene_1 |                 |          gene_2            |
	 *  ---<<<<<<<<<<----------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-----
	 *        |___________dup____________|
	 *
	 *  Then this duplication seems to disrupt gene_2:
	 *
	 *     | gene_1 |                 |          gene_2                                        |
	 *  ---<<<<<<<<<<----------------->>>><<<<<<<----------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-----
	 *        |___________dup____________||___________dup____________|
	 *
	 *  Whereas this one does not, because the duplication affects the gene
	 *  after the gene's coding region:
	 *
	 *     | gene_1 |                 |          gene_2            |
	 *  ---<<<<<<<<<<----------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-----
	 *                                     |___________dup____________|
	 *
	 *     | gene_1 |                 |          gene_2            |
	 *  ---<<<<<<<<<<----------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>--->>>>>>>>>>>>>>>>>>>>>>>>>-----
	 *                                     |___________dup____________||___________dup____________|
	 *
	 * @return true if the duplication is beyond transcript's end
	 */
	boolean beyondTranscript() {
		if (coding) {
			if (transcript.isStrandPlus()) return variant.getEnd() > transcript.getCdsEnd();
			return variant.getEnd() > transcript.getCdsStart();
		}

		return variant.getEnd() > transcript.getEnd();
	}

	@Override
	protected void effectTranscript() {
		effectNoCodon(transcript, EffectType.TRANSCRIPT_DUPLICATION);
	}

	@Override
	protected void exons() {
		if (beyondTranscript()) {
			// Is the effect of a duplication beyond transcript's end?
			// Then it probably does not have much impact
			EffectImpact impact = coding ? EffectImpact.LOW : EffectImpact.MODIFIER;
			if (exonFull > 0) effectNoCodon(transcript, EffectType.EXON_DUPLICATION, impact);
			if (exonPartial > 0) effectNoCodon(transcript, EffectType.EXON_DUPLICATION_PARTIAL, impact);
			return;
		}

		if (coding) exonsCoding();
		else exonsNoncoding();

	}

	/**
	 * One or more exons fully included (no partial overlap)
	 */
	@Override
	protected void exonsCoding() {
		codonsRefAlt();

		if (exonFull > 0) effect(transcript, EffectType.EXON_DUPLICATION, false);
		if (exonPartial > 0) effect(transcript, EffectType.EXON_DUPLICATION_PARTIAL, false);

		// Is this duplication creating a frame-shift?
		int lenDiff = cdsAlt.length() - cdsRef.length();
		if (lenDiff % 3 != 0) effect(transcript, EffectType.FRAME_SHIFT, false);
	}

	/**
	 * Effects for non-coding transcripts
	 */
	@Override
	protected void exonsNoncoding() {
		if (exonFull > 0) effectNoCodon(transcript, EffectType.EXON_DUPLICATION, EffectImpact.MODIFIER);
		if (exonPartial > 0) effectNoCodon(transcript, EffectType.EXON_DUPLICATION_PARTIAL, EffectImpact.MODIFIER);
	}

	/**
	 * Inversion does not intersect any exon
	 */
	@Override
	protected void intron() {
		effectNoCodon(transcript, EffectType.INTRON);
	}

}
