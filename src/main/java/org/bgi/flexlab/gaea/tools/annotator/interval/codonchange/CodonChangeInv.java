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
import org.bgi.flexlab.gaea.tools.annotator.interval.Marker;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant;

/**
 * Calculate codon changes produced by an inversion
 *
 * @author pcingola
 */
public class CodonChangeInv extends CodonChange {

	public CodonChangeInv(Variant variant, Transcript transcript, VariantEffects variantEffects) {
		super(variant, transcript, variantEffects);
	}

	@Override
	public void codonChange() {
		if (variant.includes(transcript)) {
			// Whole transcript inverted?
			effectNoCodon(transcript, EffectType.TRANSCRIPT_INVERSION);
		} else {
			// Part of the transcript is inverted

			// Does the inversion affect any exon?
			boolean intersectsExons = false;
			for (Exon ex : transcript) {
				if (variant.intersects(ex)) {
					intersectsExons = true;
					break;
				}
			}

			// Annotate
			if (intersectsExons) exons();
			else intron();
		}
	}

	/**
	 * One or more exons fully included (no partial overlap)
	 */
	void exons() {
		Marker cdsMarker = null;
		if (transcript.isProteinCoding()) cdsMarker = transcript.cdsMarker();

		for (Exon ex : transcript)
			if (variant.intersects(ex)) {
				EffectImpact impact = EffectImpact.LOW;

				// Is the variant affecting a coding part of the exon?
				// If so, then this is a HIGH impact effect.
				if (cdsMarker != null && variant.intersect(ex).intersects(cdsMarker)) impact = EffectImpact.HIGH;

				// Is the whole exon inverted or just part of it?
				EffectType effType = variant.includes(ex) ? EffectType.EXON_INVERSION : EffectType.EXON_INVERSION_PARTIAL;

				effectNoCodon(ex, effType, impact);
			}
	}

	/**
	 * Inversion does not intersect any exon
	 */
	void intron() {
		effectNoCodon(transcript, EffectType.INTRON);
	}

}
