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
 * Note: Splice sites donnor are defined as the last 2 bases of an intron
 * Reference: http://en.wikipedia.org/wiki/RNA_splicing
 *
 * @author pcingola
 *
 */
public class SpliceSiteAcceptor extends SpliceSite {

	private static final long serialVersionUID = -7416687954435361328L;

	public SpliceSiteAcceptor() {
		super();
		type = EffectType.SPLICE_SITE_ACCEPTOR;
	}

	public SpliceSiteAcceptor(Intron parent, int start, int end, boolean strandMinus, String id) {
		super(parent, start, end, strandMinus, id);
		type = EffectType.SPLICE_SITE_ACCEPTOR;
	}

	@Override
	public boolean intersectsCoreSpliceSite(Marker marker) {
		if (size() <= CORE_SPLICE_SITE_SIZE) return true;

		if (!getChromosomeName().equals(marker.getChromosomeName())) return false; // Not in the same chromosome? They do not intersect

		int coreStart, coreEnd;
		if (isStrandPlus()) {
			coreEnd = end;
			coreStart = coreEnd - CORE_SPLICE_SITE_SIZE + 1;
		} else {
			coreStart = start;
			coreEnd = coreStart + CORE_SPLICE_SITE_SIZE - 1;
		}

		return marker.intersects(coreStart, coreEnd);
	}
}
