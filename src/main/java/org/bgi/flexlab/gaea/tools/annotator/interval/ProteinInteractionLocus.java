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

import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect.EffectImpact;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;

import java.util.LinkedList;
import java.util.List;

/**
 * Protein interaction: An amino acid that is "in contact" with another amino acid.
 * This can be either within the same protein or interacting with another protein.
 * Evidence form PDB crystallized structures
 *
 * @author pablocingolani
 */
public abstract class ProteinInteractionLocus extends Marker {

	private static final long serialVersionUID = 1L;

	public static final boolean debug = false;

	/**
	 * Create interaction
	 */
	private static ProteinInteractionLocus factory(Transcript tr, int start, int end, Transcript trInteract, String id) {
		// Interaction type
		String geneId1 = tr.getParent().getId();
		String geneId2 = trInteract.getParent().getId();

		// Intervals may be swapped for transcript on the negative strand 
		int s = Math.min(start, end);
		int e = Math.min(start, end);

		// Same gene? => Within protein interaction
		if (geneId1.equals(geneId2)) return new ProteinStructuralInteractionLocus(tr, s, e, id);

		// Different genes? => Protein-protein interaction
		return new ProteinProteinInteractionLocus(tr, s, e, trInteract, id);

	}

	/**
	 * Create interaction. Most of the time it is only one interval, but 
	 * if introns split an amino acid, it may be more then one interval
	 */
	public static List<ProteinInteractionLocus> factory(Transcript tr, int aaPos, Transcript trInteract, String id) {
		List<ProteinInteractionLocus> list = new LinkedList<>();

		// In most cases, bases within a codon will be adjacent, but if
		// there is an intron splitting the codon, then bases will be
		// on non-contiguous positions. In such case, we need to create
		// one interaction interval for each range of contiguous bases
		// in the codon (in theory we could end up with three different
		// intervals, but that would be quite rare
		int codon2pos[] = tr.codonNumber2Pos(aaPos);

		int j = tr.isStrandPlus() ? 0 : 2;
		int start, prev, pos;
		int step = tr.isStrandPlus() ? 1 : -1;

		pos = prev = start = codon2pos[j];
		j += step;
		while (0 <= j && j <= 2) {
			pos = codon2pos[j];
			if (pos != (prev + step)) {
				// Non-contiguous, create new interval
				list.add(factory(tr, start, prev, trInteract, id));
				start = pos;
			}
			j++;
			prev = pos;
		}

		// Make sure at least one interval is created
		list.add(factory(tr, start, pos, trInteract, id));
		return list;
	}

	public ProteinInteractionLocus() {
		super();
	}

	public ProteinInteractionLocus(Transcript parent, int start, int end, String id) {
		super(parent, start, end, false, id);

	}

	@Override
	public ProteinInteractionLocus cloneShallow() {
		ProteinInteractionLocus clone = (ProteinInteractionLocus) super.cloneShallow();
		return clone;
	}

	/**
	 * Calculate the effect of this variant
	 */
	@Override
	public boolean variantEffect(Variant variant, VariantEffects variantEffects) {
		if (!intersects(variant)) return false;// Sanity check
		variantEffects.add(variant, this, type, EffectImpact.HIGH, "");
		return true;
	}

}
