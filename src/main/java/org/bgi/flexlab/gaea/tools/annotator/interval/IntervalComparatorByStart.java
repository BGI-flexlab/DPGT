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

import java.util.Comparator;

/**
 * Compare intervals by start position
 * @author pcingola
 *
 */
public class IntervalComparatorByStart implements Comparator<Marker> {

	int order = 1;

	public IntervalComparatorByStart() {
		super();
	}

	public IntervalComparatorByStart(boolean reverse) {
		super();
		if (reverse) order = -1;
	}

	@Override
	public int compare(Marker i1, Marker i2) {
		// Compare chromosome
		if ((i1.getChromosomeNum() == 0) || (i2.getChromosomeNum() == 0)) { // Use string version?
			// Chromosome by string
			int c = i1.getChromosomeName().compareTo(i2.getChromosomeName());
			if (c != 0) return order * c;
		} else {
			// Use numeric version
			if (i1.getChromosomeNum() > i2.getChromosomeNum()) return order;
			if (i1.getChromosomeNum() < i2.getChromosomeNum()) return -order;
		}

		// Start
		if (i1.getStart() > i2.getStart()) return order;
		if (i1.getStart() < i2.getStart()) return -order;

		// End
		if (i1.getEnd() > i2.getEnd()) return order;
		if (i1.getEnd() < i2.getEnd()) return -order;

		// Compare by ID
		if ((i1.getId() == null) && (i2.getId() == null)) return 0;
		if ((i1.getId() != null) && (i2.getId() == null)) return -1;
		if ((i1.getId() == null) && (i2.getId() != null)) return 1;
		return i1.getId().compareTo(i2.getId());
	}
}
