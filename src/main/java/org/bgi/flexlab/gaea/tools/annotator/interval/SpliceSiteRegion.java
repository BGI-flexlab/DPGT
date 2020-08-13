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
 * Interval for a splice site acceptor
 *
 * From Sequence Ontology: A sequence variant in which a change has occurred
 * within the region of the splice site, either within 1-3 bases of the exon
 * or 3-8 bases of the intron.
 *
 * @author pcingola
 */
public class SpliceSiteRegion extends SpliceSite {

	private static final long serialVersionUID = -7416687954435361328L;

	public SpliceSiteRegion() {
		super();
		type = EffectType.SPLICE_SITE_REGION;
	}

	public SpliceSiteRegion(Exon parent, int start, int end, boolean strandMinus, String id) {
		super(parent, start, end, strandMinus, id);
		type = EffectType.SPLICE_SITE_REGION;
	}

	public SpliceSiteRegion(Intron parent, int start, int end, boolean strandMinus, String id) {
		super(parent, start, end, strandMinus, id);
		type = EffectType.SPLICE_SITE_REGION;
	}

	@Override
	public boolean intersectsCoreSpliceSite(Marker marker) {
		return false;
	}

	public boolean isExonPart() {
		return parent instanceof Exon;
	}

	public boolean isIntronPart() {
		return parent instanceof Intron;
	}
}
