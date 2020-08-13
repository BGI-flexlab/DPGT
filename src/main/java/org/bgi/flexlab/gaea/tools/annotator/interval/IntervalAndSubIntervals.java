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

import java.util.*;

/**
 * Interval that contains sub intervals.
 *
 * @author pcingola
 *
 */
public class IntervalAndSubIntervals<T extends Marker> extends Marker implements Iterable<T> {

	private static final long serialVersionUID = 1636197649250882952L;
	HashMap<String, T> subIntervals;
	ArrayList<T> sorted;
	ArrayList<T> sortedStrand;

	public IntervalAndSubIntervals() {
		super();
		reset();
	}

	public IntervalAndSubIntervals(Marker parent, int start, int end, boolean strandMinus, String id) {
		super(parent, start, end, strandMinus, id);
		reset();
	}

	/**
	 * Add a subinterval
	 */
	public synchronized void add(T t) {
		if (subIntervals.put(t.getId(), t) != null) {
			// Keys should be unique
			throw new RuntimeException(t.getClass().getSimpleName() //
					+ " '" + t.getId() + "' is already in " //
					+ this.getClass().getSimpleName() //
					+ " '" + id + "'" //
			);
		}

		// Sort is no longer valid
		invalidateSorted();
	}

	/**
	 * Add all intervals
	 */
	public void addAll(Iterable<T> ts) {
		for (T t : ts)
			add(t);
		invalidateSorted();
	}

	/**
	 * Add all markers
	 */
	@SuppressWarnings("unchecked")
	public void addAll(Markers markers) {
		for (Marker m : markers)
			add((T) m);
		invalidateSorted();
	}

	/**
	 * Apply a variant.
	 */
	@Override
	@SuppressWarnings("unchecked")
	public IntervalAndSubIntervals<T> apply(Variant variant) {
		if (!shouldApply(variant)) return this;

		// Apply to Marker
		IntervalAndSubIntervals<T> newMarker = (IntervalAndSubIntervals<T>) super.apply(variant);
		if (newMarker == null) return null;

		// Now apply to all sub-markers
		newMarker.reset();
		for (T m : this) {
			// Whole marker duplication
			if (variant.isDup() && variant.includes(m)) {
				// We need to add another copy of the marker: One without any modification and one with coordinates updated
				T mcopy = (T) m.cloneShallow();
				mcopy.setId(mcopy.getId() + ".dup"); // Change ID
				mcopy.setParent(newMarker);
				newMarker.add(mcopy);
			}

			// Apply variant to marker
			T mcopy = (T) m.apply(variant);

			// Do not add if interval is completely removed
			if (mcopy != null) {
				// Make sure we don't modify the original subinterval.
				if (mcopy == m) mcopy = (T) m.cloneShallow();

				mcopy.setParent(newMarker);
				newMarker.add(mcopy);
			}
		}

		return newMarker;
	}

	@SuppressWarnings("unchecked")
	@Override
	public IntervalAndSubIntervals<T> clone() {
		IntervalAndSubIntervals<T> copy = (IntervalAndSubIntervals<T>) super.clone();
		copy.reset();

		for (T m : this) {
			T mcopy = (T) m.clone();
			mcopy.setParent(copy);
			copy.add(mcopy);
		}

		return copy;
	}

	@SuppressWarnings("unchecked")
	@Override
	public IntervalAndSubIntervals<T> cloneShallow() {
		IntervalAndSubIntervals<T> clone = (IntervalAndSubIntervals<T>) super.cloneShallow();
		clone.reset();
		return clone;
	}

	/**
	 * Is 'id' in the subintervals?
	 */
	public boolean containsId(String id) {
		return subIntervals.containsKey(id);
	}

	/**
	 * Obtain a subinterval
	 */
	public T get(String id) {
		return subIntervals.get(id);
	}

	/**
	 * Invalidate sorted collections
	 */
	protected void invalidateSorted() {
		sorted = sortedStrand = null;
	}

	@Override
	public Iterator<T> iterator() {
		return subIntervals().iterator();
	}

	/**
	 * A list of all markers in this transcript
	 */
	public Markers markers() {
		Markers markers = new Markers();
		markers.addAll(subIntervals());
		return markers;
	}

	public int numChilds() {
		return (subIntervals != null ? subIntervals.size() : 0);
	}

	/**
	 * Query all genomic regions that intersect 'marker'
	 */
	@Override
	public Markers query(Marker marker) {
		Markers markers = new Markers();

		for (Marker m : this) {
			if (m.intersects(marker)) {
				markers.add(m);

				Markers subMarkers = m.query(marker);
				if (subMarkers != null) markers.add(subMarkers);
			}
		}

		return markers;
	}

	/**
	 * Remove a subinterval
	 */
	public synchronized void remove(T t) {
		subIntervals.remove(t.getId());
		invalidateSorted();
	}

	/**
	 * Remove all intervals
	 */
	public synchronized void reset() {
		subIntervals = new HashMap<String, T>();
		invalidateSorted();
	}

	@Override
	public void setStrandMinus(boolean strandMinus) {
		this.strandMinus = strandMinus;

		// Change all subintervals
		for (T t : this)
			t.setStrandMinus(strandMinus);

		invalidateSorted(); // These are no longer correct
	}

	@Override
	public void shiftCoordinates(int shift) {
		super.shiftCoordinates(shift);

		for (T m : subIntervals())
			m.shiftCoordinates(shift);
	}

	/**
	 * Return a collection of sub intervals sorted by natural order
	 */
	public synchronized List<T> sorted() {
		if (sorted != null) return sorted;
		ArrayList<T> sorted = new ArrayList<T>();
		sorted.addAll(subIntervals());
		Collections.sort(sorted);

		this.sorted = sorted;
		return sorted;
	}

	/**
	 * Return a collection of sub intervals sorted by start position (if strand is >= 0) or
	 * by reverse end position (if strans < 0)
	 */
	public synchronized List<T> sortedStrand() {
		if (sortedStrand != null) return sortedStrand;

		ArrayList<T> sortedStrand = new ArrayList<T>();
		sortedStrand.addAll(subIntervals());

		if (isStrandPlus()) Collections.sort(sortedStrand, new IntervalComparatorByStart()); // Sort by start position
		else Collections.sort(sortedStrand, new IntervalComparatorByEnd(true)); // Sort by end position (reversed)

		this.sortedStrand = sortedStrand;
		return sortedStrand;
	}

	/**
	 * Return a collection of sub intervals
	 */
	public Collection<T> subIntervals() {
		return subIntervals.values();
	}

}
