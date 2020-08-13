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


/**
 * Interval for in intergenic region
 *
 * @author pcingola
 *
 */
public class Intergenic extends Marker {

	private static final long serialVersionUID = -2487664381262354896L;

	String name;

	public Intergenic() {
		super();
		type = EffectType.INTERGENIC;
		name = "";
	}

	public Intergenic(Chromosome parent, int start, int end, boolean strandMinus, String id, String name) {
		super(parent, start, end, strandMinus, id);
		type = EffectType.INTERGENIC;
		this.name = name;
	}

	@Override
	public Intergenic cloneShallow() {
		Intergenic clone = (Intergenic) super.cloneShallow();
		clone.name = name;
		return clone;
	}

	public String getName() {
		return name;
	}

}
