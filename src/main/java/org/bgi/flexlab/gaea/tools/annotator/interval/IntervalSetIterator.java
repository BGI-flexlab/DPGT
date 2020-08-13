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

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Iterate over intervals.
 * Note: The returned set is recycled on every iteration
 *  
 * @author pcingola
 *
 */
public class IntervalSetIterator implements Iterator<Set<Marker>>, Iterable<Set<Marker>> {

	boolean firstIteration = true;
	boolean newChromosome = true;
	HashSet<Marker> openIntervals; // Set of intervals in out iterator
	HashSet<Marker> openIntervalsOut; // Same as 'openIntervals', nut this one is exposed to the "outside world" (i.e. return in 'next()' function)
	Markers intervalsByStart, intervalsByEnd; // Sort intervals by start and end position
	Iterator<Marker> itStart, itEnd; // Iterate on intervals sorted by start/end position
	Marker intStart, intEnd; // Current intervals in iteration (sorted by start and end)
	int start, end;
	int i;
	String chromo; // Current chromosome

	public IntervalSetIterator(Markers intervals) {
		openIntervals = new HashSet<Marker>();
		openIntervalsOut = new HashSet<Marker>();

		intervalsByStart = new Markers();
		intervalsByStart.add(intervals);
		intervalsByStart.sort(false, false);

		intervalsByEnd = new Markers();
		intervalsByEnd.add(intervals);
		intervalsByEnd.sort(true, false);

		itStart = intervalsByStart.iterator();
		itEnd = intervalsByEnd.iterator();
		i = start = end = -1;
	}

	/**
	 * Is 'chr' equal to current chromosome?
	 * @param chr
	 * @return
	 */
	boolean equalsChromo(String chr) {
		if( (chr == null) && (chromo == null) ) return true;
		if( (chr == null) || (chromo == null) ) return false;
		return chr.equals(chromo);
	}

	@Override
	public boolean hasNext() {
		return itStart.hasNext() || itEnd.hasNext() || (openIntervals.size() > 0);
	}

	/**
	 * Iteration
	 * This is not a trivial process (or may be I overly complicated it :-)
	 */
	void iterate() {
		if( firstIteration ) { // First iteration: Initialize
			iterateInit();
		} else if( newChromosome ) { // Finished iterating on a chromosome and jumping to the next? => we need to do some special stuff
			if( equalsChromo(intEnd.getChromosomeName()) ) iterateNext(); // Still some intervals from the old chromosome in the 'End' list => finish iterating on them
			else iterateInitNewChromo(); // Initialize for new chromosome
		} else if( itStart.hasNext() || itEnd.hasNext() ) { // Has more items to iterate?
			iterateNext(); // Prepare next
		} else if( openIntervals.size() > 0 ) { // No more items to iterate, but some remain in the 'openIntervals'? => remove them one by one
			iterateFinish();
		}
	}

	/**
	 * Finish iteration steps
	 */
	void iterateFinish() {
		Interval interv = openIntervals.iterator().next();
		openIntervals.remove(interv);

	}

	/**
	 * Initialize iteration
	 */
	void iterateInit() {
		intStart = itStart.next();
		openIntervals.add(intStart);
		start = intStart.getStart();
		chromo = intStart.getChromosomeName();

		intEnd = itEnd.next();
		end = intEnd.getEnd();
		openIntervals.add(intEnd);

		i = intStart.getStart();
		firstIteration = false;
		newChromosome = false;
	}

	/**
	 * Initialize iteration when jumping to a new chromosome
	 */
	void iterateInitNewChromo() {
		chromo = intStart.getChromosomeName();
		i = start = end = -1;
		intEnd = null;

		openIntervals.add(intStart);
		start = intStart.getStart();

		//intEnd = itEnd.next();
		end = intEnd.getEnd();
		openIntervals.add(intEnd);

		i = Math.min(start, end);
		newChromosome = false;
	}

	/**
	 * Next iteration step
	 */
	void iterateNext() {
		// Need to iterate on intervalsByStart?
		if( (intStart != null) && (intStart.getStart() <= i) ) {

			if( newChromosome ) {
				// Nothing to do
			} else if( itStart.hasNext() ) {
				intStart = itStart.next();

				// Only proceed if it's in the same chromosome
				if( equalsChromo(intStart.getChromosomeName()) ) {
					openIntervals.add(intStart);
					start = intStart.getStart();
					i = Math.min(start, end);
					return;
				} else newChromosome = true;
			} else {
				intStart = null;
				start = Integer.MAX_VALUE;
			}
		}

		// Need to iterate on intervalsByStart?
		if( ((intEnd != null) && (intEnd.getEnd() <= i)) || newChromosome || (intStart == null) ) {
			if( itEnd.hasNext() ) {
				openIntervals.remove(intEnd);
				intEnd = itEnd.next();
				if( equalsChromo(intEnd.getChromosomeName()) ) {
					end = intEnd.getEnd();
					i = Math.min(start, end);
				}
				return;
			} else {
				openIntervals.remove(intEnd);
				intEnd = null;
				return;
			}
		}

	}

	@Override
	public Iterator<Set<Marker>> iterator() {
		return this;
	}

	@Override
	public Set<Marker> next() {
		if( !hasNext() ) return null;
		iterate();
		return openIntervals;
	}

	@Override
	public void remove() {
		throw new RuntimeException("Unimplmented method!");
	}

	@Override
	public String toString() {
		return "i:" + i + "\tintStart: " + intStart + "\tintEnd: " + intEnd + "\topenIntervals: " + openIntervals;
	}
}
