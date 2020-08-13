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

/**
 * A variant respect to non-reference (e.g. comparing cancer vs. somatic tissue).
 *
 * @author pcingola
 */
public class VariantNonRef extends Variant {

	private static final long serialVersionUID = 1L;
	Variant variantRef;

	public VariantNonRef() {
		super();
	}

	public VariantNonRef(Variant variant, Variant variantRef) {
		super(variant.getParent(), variant.getStart(), variant.getReference(), variant.getAlt(), variant.getId());
		genotype = variant.getGenotype();
		this.variantRef = variantRef;
	}

	@Override
	public String getGenotype() {
		return genotype + "-" + variantRef.getGenotype();
	}

	public Variant getVariantRef() {
		return variantRef;
	}

	@Override
	public boolean isNonRef() {
		return true;
	}

	@Override
	public Variant realignLeft() {
		// Realigning in cancer samples is not trivial: What happens if one realigns and the other doesn't?
		// For now, do not realign
		return this;
	}

	public void setVariantRef(Variant variantRef) {
		this.variantRef = variantRef;
	}

	@Override
	public String toString() {
		String valt = super.toString();
		String vref = variantRef.toString();
		return valt + "-" + vref;
	}
}
