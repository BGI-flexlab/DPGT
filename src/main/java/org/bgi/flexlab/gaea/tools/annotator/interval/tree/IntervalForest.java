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

import org.bgi.flexlab.gaea.tools.annotator.interval.Chromosome;
import org.bgi.flexlab.gaea.tools.annotator.interval.Marker;
import org.bgi.flexlab.gaea.tools.annotator.interval.Markers;
import org.bgi.flexlab.gaea.tools.annotator.util.Gpr;

import java.io.Serializable;
import java.util.*;

/**
 * A set of interval trees (e.g. one per chromosome, one per transcript ID, etc)
 *
 * @author pcingola
 */
public class IntervalForest implements Serializable, Iterable<Itree> {

	private static final long serialVersionUID = 1L;

	boolean debug;
	HashMap<String, Itree> forest;

	public IntervalForest() {
		forest = new HashMap<String, Itree>();
	}

	public IntervalForest(Markers markers) {
		forest = new HashMap<String, Itree>();
		add(markers);
	}

	/**
	 * Add all intervals
	 */
	public void add(Collection<? extends Marker> intervals) {
		for (Marker i : intervals)
			add(i);
	}

	/**
	 * Add an interval
	 */
	public void add(Marker interval) {
		if (interval == null) return;
		String chName = Chromosome.simpleName(interval.getChromosomeName());
		getOrCreateTreeChromo(chName).add(interval); // Add interval to tree
	}

	/**
	 * Add all intervals
	 */
	public void add(Markers intervals) {
		for (Marker i : intervals)
			add(i);
	}

	/**
	 * Build all trees
	 */
	public void build() {
		for (String key : forest.keySet()) {
			if (debug) Gpr.debug("Building interval tree for '" + key + "'");
			Itree tree = forest.get(key);
			tree.build();
		}
	}

	/**
	 * Get (or create) an interval tree for ID
	 */
	public Itree getOrCreateTree(String id) {
		// Retrieve (or create) interval tree
		Itree itree = forest.get(id);
		if (itree == null) {
			itree = newItree();
			itree.build();
			forest.put(id, itree);
		}

		return itree;
	}

	/**
	 * Get (or create) an interval tree based for "chromo" (chromosome name)
	 */
	public Itree getOrCreateTreeChromo(String chromo) {
		return getOrCreateTree(Chromosome.simpleName(chromo));
	}

	/**
	 * Get an interval tree using an ID
	 */
	public Itree getTree(String key) {
		return forest.get(key);
	}

	/**
	 * Get an interval tree using a chromosome name
	 */
	public Itree getTreeChromo(String chromo) {
		return forest.get(Chromosome.simpleName(chromo));
	}

	/**
	 * Is the tree 'chromo' available?
	 */
	public boolean hasTree(String chromo) {
		return getTreeChromo(chromo) != null;
	}

	/**
	 * Return the intersection of 'markers' and this IntervalForest
	 *
	 * For each marker 'm' in 'markers'
	 * 		- query the tree to get all markers intersecting 'm'
	 * 		- create a new interval which is the intersection of 'm' with all the resutls from the previous query.
	 */
	public Markers intersect(Markers markers) {
		Markers result = new Markers();

		// Add all intersecting intervals
		for (Marker mm : markers) {
			Markers query = query(mm);
			if (query != null) {
				for (Marker mq : query) {
					// Intersection between 'mm' and 'mq'
					int start = Math.max(mq.getStart(), mm.getStart());
					int end = Math.max(mq.getEnd(), mm.getEnd());
					Marker mintq = new Marker(mq.getParent(), start, end, mq.isStrandMinus(), "");

					// Add intersection result
					result.add(mintq);
				}
			}
		}

		return result;
	}

	@Override
	public Iterator<Itree> iterator() {
		return forest.values().iterator();
	}

	public Collection<String> keySet() {
		return forest.keySet();
	}

	/**
	 * Create new tree.
	 * In oder to change the implementation, only this method should be changed.
	 */
	protected Itree newItree() {
		return new IntervalTree();
	}

	/**
	 * Query all intervals that intersect with 'interval'
	 */
	public Markers query(Marker marker) {
		return getOrCreateTreeChromo(marker.getChromosomeName()).query(marker);
	}

	/**
	 * Query all intervals that intersect with any interval in 'intervals'
	 */
	public Markers query(Markers marker) {
		Markers ints = new Markers();

		// Add all intersecting intervals
		for (Marker i : marker)
			ints.add(query(i));

		return ints;
	}

	/**
	 * Query unique intervals that intersect with any interval in 'markers'
	 * I.e.: Return a set of intervals that intersects (at least once) with any interval in 'markers'
	 */
	public Markers queryUnique(Markers markers) {
		HashSet<Marker> uniqueMarkers = new HashSet<Marker>();

		// Add all intersecting intervals
		for (Marker q : markers) {
			Markers results = query(q); // Query

			for (Marker r : results)
				// Add all results
				uniqueMarkers.add(r);
		}

		// Create markers
		Markers ints = new Markers();
		for (Marker r : uniqueMarkers)
			ints.add(r);

		return ints;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public int size() {
		int size = 0;
		for (Itree it : forest.values())
			size += it.size();
		return size;
	}

	/**
	 * Obtain all intervals that intersect with 'marker.start'
	 */
	public Markers stab(Marker marker) {
		return stab(marker.getChromosomeName(), marker.getStart());
	}

	/**
	 * Obtain all intervals that intersect with 'point'
	 */
	public Markers stab(String chromo, int point) {
		return getOrCreateTreeChromo(chromo).stab(point);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		ArrayList<String> keys = new ArrayList<>();
		keys.addAll(forest.keySet());
		Collections.sort(keys);

		for (String key : keys) {
			Itree tree = getOrCreateTreeChromo(key);
			sb.append(key + "\tsize:" + tree.size() + "\tin_sync: " + tree.isInSync() + "\n");
		}

		return sb.toString();
	}

}
