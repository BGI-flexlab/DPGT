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

import org.bgi.flexlab.gaea.tools.annotator.interval.tree.IntervalForest;

/**
 * A translocation consisting of two endpoints
 *
 * @author pcingola
 */
public class VariantTranslocation extends Variant {

	private static final long serialVersionUID = 1L;
	Marker endPoint;
	boolean left; // Is endPoint oriented to the left?
	boolean before; // Is endPoint before ALT?

	public VariantTranslocation() {
		super();
	}

	public VariantTranslocation(Marker parent, int start, String ref, String alt, Chromosome chrTr, int startTr, boolean left, boolean before) {
		super(parent, start, ref, alt);
		endPoint = new Marker(chrTr, startTr, startTr);
		endPoint.setStrandMinus(left);
		this.left = left;
		this.before = before;
		variantType = VariantType.BND;
	}

	public Marker getEndPoint() {
		return endPoint;
	}

	public boolean isBefore() {
		return before;
	}

	@Override
	public boolean isBnd() {
		return true;
	}

	public boolean isLeft() {
		return left;
	}

	@Override
	public boolean isStructural() {
		return true;
	}

	/**
	 * Return a collection of intervals that intersect both ends of this variant
	 */
	@Override
	public Markers query(IntervalForest intervalForest) {
		Markers res1 = intervalForest.query(this);
		Markers res2 = intervalForest.query(endPoint);
		res1.add(res2);
		return res1;
	}

	@Override
	public Variant realignLeft() {
		// Do not realign translocations
		return this;
	}

	@Override
	public String toString() {
		String sep = left ? "]" : "[";
		String trPos = sep + endPoint.getChromosomeName() + ":" + endPoint.getStart() + sep;

		return "chr" + getChromosomeName() //
				+ ":" + start //
				+ "_" + getReference() //
				+ "/" //
				+ (before ? trPos + getAlt() : getAlt() + trPos) //
				;
	}

}
