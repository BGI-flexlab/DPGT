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

import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;

/**
 * Interval for a splice site
 *
 * Reference: http://en.wikipedia.org/wiki/RNA_splicing
 *
 * Spliceosomal introns often reside in eukaryotic protein-coding genes. Within the
 * intron, a 3' splice site, 5' splice site, and branch site are required for splicing.
 * The 5' splice site or splice donor site includes an almost invariant sequence GU
 * at the 5' end of the intron, within a larger, less highly conserved consensus region.
 * The 3' splice site or splice acceptor site terminates the intron with an almost
 * invariant AG sequence. Upstream (5'-ward) from the AG there is a region high in
 * pyrimidines (C and U), or polypyrimidine tract. Upstream from the polypyrimidine
 * tract is the branch point, which includes an adenine nucleotide.
 *
 * @author pcingola
 *
 */
public abstract class SpliceSite extends Marker {

	private static final long serialVersionUID = 1636197649250882952L;

	public static final int CORE_SPLICE_SITE_SIZE = 2;
	public static final int SPLICE_REGION_EXON_SIZE = 3;
	public static final int SPLICE_REGION_INTRON_MIN = 3;
	public static final int SPLICE_REGION_INTRON_MAX = 8;

	public SpliceSite() {
		super();
	}

	public SpliceSite(Exon parent, int start, int end, boolean strandMinus, String id) {
		super(parent, start, end, strandMinus, id);
	}

	public SpliceSite(Intron parent, int start, int end, boolean strandMinus, String id) {
		super(parent, start, end, strandMinus, id);
	}

	/**
	 * Core splice sites are defined as CORE_SPLICE_SITE_SIZE bases after exon end
	 * or before exon begins. Usually CORE_SPLICE_SITE_SIZE is 2 bases.
	 * Other spice sites are considered "non-core".
	 */
	public abstract boolean intersectsCoreSpliceSite(Marker marker);

	/**
	 * Splice sites are no included in Exons, by definition.
	 */
	@Override
	protected boolean isShowWarningIfParentDoesNotInclude() {
		return false;
	}

	@Override
	public boolean variantEffect(Variant variant, VariantEffects variantEffects) {
		if (!intersects(variant)) return false; // Sanity check
		variantEffects.addEffectType(variant, this, type);
		return true;
	}
}
