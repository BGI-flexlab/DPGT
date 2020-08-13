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
 * Regulatory elements
 *
 * @author pablocingolani
 */
public class Regulation extends Marker {

	private static final long serialVersionUID = -5607588295343642199L;

	String cellType = "";
	String name = "";

	public Regulation() {
		super();
		type = EffectType.REGULATION;
	}

	public Regulation(Marker parent, int start, int end, boolean strandMinus, String id, String name, String cellType) {
		super(parent, start, end, strandMinus, id);
		type = EffectType.REGULATION;
		this.name = name;
		this.cellType = cellType;
	}

	@Override
	public Regulation cloneShallow() {
		Regulation clone = (Regulation) super.cloneShallow();
		clone.cellType = cellType;
		clone.name = name;
		return clone;
	}

	public String getCellType() {
		return cellType;
	}

	public String getName() {
		return name;
	}

	@Override
	public String toString() {
		return getChromosomeName() + "\t" + start + "-" + end //
				+ " " //
				+ type + ((name != null) && (!name.isEmpty()) ? " '" + name + "'" : "");
	}

	/**
	 * Calculate the effect of this seqChange
	 */
	@Override
	public boolean variantEffect(Variant variant, VariantEffects variantEffects) {
		if (!intersects(variant)) return false; // Sanity check
		EffectType effType = EffectType.REGULATION;
		variantEffects.add(variant, this, effType, "");
		return true;
	}

}
