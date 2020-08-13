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
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant;
import org.bgi.flexlab.gaea.tools.annotator.util.Gpr;

import java.util.List;

/**
 * Calculate codon changes produced by a 'mixed' variant
 *
 * Essentially every 'mixed' variant can be represented as a concatenation of a SNP/MNP + an INS/DEL
 *
 * @author pcingola
 */
public class CodonChangeMixed extends CodonChangeMnp {

	public static boolean debug = false;

	int oldCodonCdsStart = -1;
	int oldCodonCdsEnd = -1;

	Variant mnp;
	Variant indel;
	CodonChange codonChangeMnp;
	CodonChange codonChangeIndel;
	VariantEffects variantEffectsOri;

	public CodonChangeMixed(Variant variant, Transcript transcript, VariantEffects variantEffects) {
		super(variant, transcript, variantEffects);
		returnNow = false;
		requireNetCdsChange = false;

		// Decompose mixed variant into 'SNP/MNP' + 'INS/DEL'
		int minLen = Math.min(variant.getReference().length(), variant.getAlt().length());

		String ref = variant.getReference();
		String refMnp = ref.substring(0, minLen);
		String refIndel = ref.substring(minLen);

		String alt = variant.getAlt();
		String altMnp = alt.substring(0, minLen);
		String altIndel = alt.substring(minLen);

		mnp = new Variant(variant.getChromosome(), variant.getStart(), refMnp, altMnp, variant.getId());
		indel = new Variant(variant.getChromosome(), variant.getStart() + minLen, refIndel, altIndel, variant.getId());

		// Create codon changes
		variantEffectsOri = variantEffects;
		this.variantEffects = new VariantEffects();
		codonChangeMnp = CodonChange.factory(mnp, transcript, this.variantEffects);
		codonChangeIndel = CodonChange.factory(indel, transcript, this.variantEffects);
	}

	@Override
	public void codonChange() {
		codonOldNew();

		codonChangeMnp.codonChange();
		codonChangeIndel.codonChange();

		// Set highest impact variant effect
		if (variantEffects.isEmpty()) return; // Nothing to do

		variantEffects.sort();
		VariantEffect varEff = variantEffects.get(0);

		if (debug) {
			Gpr.debug("Mixed variant:" + variant + "\n\t\tSNP/MNP : " + mnp + "\n\t\tInDel   : " + indel + "\n\t\tEffects : ");
			for (VariantEffect ve : variantEffects)
				System.err.println("\t\t\t" + ve.toStringSimple(true));
		}

		// Add main effect
		varEff = effect(varEff.getMarker(), varEff.getEffectType(), false);

		// Add 'additional' effects
		for (int i = 0; i < variantEffects.size(); i++) {
			List<EffectType> effTypes = variantEffects.get(i).getEffectTypes();
			for (int j = 0; j < effTypes.size(); j++) {
				EffectType effType = effTypes.get(j);
				if (!varEff.hasEffectType(effType)) varEff.addEffect(effType);
			}
		}

		variantEffectsOri.add(varEff);
	}

	void codonNum() {
		if (transcript.isStrandPlus()) {
			codonStartNum = codonChangeMnp.codonStartNum;
			codonStartIndex = codonChangeMnp.codonStartIndex;
		} else {
			codonStartNum = codonChangeIndel.codonStartNum;
			codonStartIndex = codonChangeIndel.codonStartIndex;
		}
	}

}
