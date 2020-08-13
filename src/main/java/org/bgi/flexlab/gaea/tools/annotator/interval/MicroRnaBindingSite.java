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
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;

/**
 * miRna binding site (usually this was predicted by some algorithm)
 *
 * @author pcingola
 */
public class MicroRnaBindingSite extends Marker {

	private static final long serialVersionUID = -9089500641817245554L;

	double pValue;

	public MicroRnaBindingSite() {
		super();
		type = EffectType.MICRO_RNA;
	}

	public MicroRnaBindingSite(Marker parent, int start, int end, boolean strandMinus, String id, double pValue) {
		super(parent, start, end, strandMinus, id);
		this.pValue = pValue;
		type = EffectType.MICRO_RNA;
	}

	@Override
	public MicroRnaBindingSite cloneShallow() {
		MicroRnaBindingSite clone = (MicroRnaBindingSite) super.cloneShallow();
		clone.pValue = pValue;
		return clone;
	}

	@Override
	public boolean variantEffect(Variant variant, VariantEffects changeEffects) {
		if (!intersects(variant)) return false; // Sanity check
		changeEffects.add(variant, this, EffectType.MICRO_RNA, "" + pValue);
		return true;
	}

}
