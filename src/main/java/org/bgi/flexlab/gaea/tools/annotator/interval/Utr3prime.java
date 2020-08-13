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

import org.bgi.flexlab.gaea.tools.annotator.effect.EffectType;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant.VariantType;

/**
 * Interval for a UTR (5 prime UTR and 3 prime UTR
 *
 * @author pcingola
 *
 */
public class Utr3prime extends Utr {

	private static final long serialVersionUID = 5688641008301281991L;

	public Utr3prime() {
		super();
		type = EffectType.UTR_3_PRIME;
	}

	public Utr3prime(Exon parent, int start, int end, boolean strandMinus, String id) {
		super(parent, start, end, strandMinus, id);
		type = EffectType.UTR_3_PRIME;
	}

	@Override
	public boolean isUtr3prime() {
		return true;
	}

	@Override
	public boolean isUtr5prime() {
		return false;
	}

	/**
	 * Calculate distance from beginning of 3'UTRs
	 */
	@Override
	int utrDistance(Variant variant, Transcript tr) {
		int cdsEnd = tr.getCdsEnd();
		if (cdsEnd < 0) return -1;

		if (isStrandPlus()) return variant.getStart() - cdsEnd;
		return cdsEnd - variant.getEnd();
	}

	@Override
	public boolean variantEffect(Variant variant, VariantEffects variantEffects) {
		if (!intersects(variant)) return false;

		if (variant.includes(this) && (variant.getVariantType() == VariantType.DEL)) {
			variantEffects.addEffectType(variant, this, EffectType.UTR_3_DELETED); // A UTR was removed entirely
			return true;
		}

		Transcript tr = (Transcript) findParent(Transcript.class);
		int distance = utrDistance(variant, tr);

		VariantEffect variantEffect = new VariantEffect(variant);
		variantEffect.set(this, type, type.effectImpact(), distance >= 0 ? distance + " bases from CDS" : "");
		variantEffect.setDistance(distance);
		variantEffects.add(variantEffect);

		return true;
	}
}
