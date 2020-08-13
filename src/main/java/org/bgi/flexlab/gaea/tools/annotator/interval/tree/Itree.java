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
package org.bgi.flexlab.gaea.tools.annotator.interval.tree;

import org.bgi.flexlab.gaea.tools.annotator.interval.Interval;
import org.bgi.flexlab.gaea.tools.annotator.interval.Marker;
import org.bgi.flexlab.gaea.tools.annotator.interval.Markers;

/**
 * Interval tree interface
 */
public interface Itree extends Iterable<Marker> {

	/**
	 * Add an interval object to the interval tree's list
	 */
	public void add(Marker interval);

	/**
	 * Add all intervals to interval tree's list
	 */
	public void add(Markers markers);

	/**
	 * Build the interval tree to reflect the list of intervals.
	 * Must not run if this is currently in sync
	 */
	public void build();

	public Markers getIntervals();

	public boolean isEmpty();

	/**
	 * Is the tree 'in sync'?
	 * If false, the tree must be 'build()' before the next query
	 */
	public boolean isInSync();

	/**
	 * Load intervals from file
	 */
//	public void load(String fileName, Genome genome);

	/**
	 * Perform an interval query, returning the intervals that
	 * intersect with 'interval'
	 *
	 * @return All intervals that intersect 'interval'
	 */
	public Markers query(Interval interval);

	/**
	 * Size: number of entries in this tree
	 */
	public int size();

	/**
	 * Perform a stabbing query, returning the interval objects
	 * @return All intervals intersecting 'point'
	 */
	public Markers stab(int point);

}
