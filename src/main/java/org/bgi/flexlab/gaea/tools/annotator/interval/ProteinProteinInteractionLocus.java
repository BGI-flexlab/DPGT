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
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect.EffectImpact;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;

/**
 * Protein interaction: An amino acid that is "in contact" with another amino acid
 * within the same protein. Evidence form PDB crystallized structures
 *
 * @author pablocingolani
 */
public class ProteinProteinInteractionLocus extends ProteinInteractionLocus {

	private static final long serialVersionUID = -2111056845758378137L;

	Transcript trInteract;

	public ProteinProteinInteractionLocus() {
		super();
		type = EffectType.PROTEIN_PROTEIN_INTERACTION_LOCUS;
	}

	public ProteinProteinInteractionLocus(Transcript parent, int start, int end, Transcript trInteract, String id) {
		super(parent, start, end, id);
		this.trInteract = trInteract;
		type = EffectType.PROTEIN_PROTEIN_INTERACTION_LOCUS;
	}

	@Override
	public ProteinProteinInteractionLocus cloneShallow() {
		ProteinProteinInteractionLocus clone = (ProteinProteinInteractionLocus) super.cloneShallow();
		clone.trInteract = trInteract;
		return clone;
	}

	/**
	 * Calculate the effect of this variant
	 */
	@Override
	public boolean variantEffect(Variant variant, VariantEffects variantEffects) {
		if (!intersects(variant)) return false;// Sanity check
		variantEffects.add(variant, this, EffectType.PROTEIN_PROTEIN_INTERACTION_LOCUS, EffectImpact.HIGH, "");
		return true;
	}

}
