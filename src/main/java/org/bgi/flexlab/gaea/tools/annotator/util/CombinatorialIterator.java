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
package org.bgi.flexlab.gaea.tools.annotator.util;

import java.util.Iterator;

/**
 * Generate all possible 'count' combinations
 *
 * @author pcingola
 */
public class CombinatorialIterator implements Iterable<int[]>, Iterator<int[]> {

	int next[];
	int max[];
	int min[];
	boolean inc, finished;

	public CombinatorialIterator(int size) {
		next = new int[size];
		max = new int[size];
		min = new int[size];
	}

	@Override
	public boolean hasNext() {
		if (inc) inc();
		return !finished;
	}

	public void inc() {
		inc = false;

		for (int i = 0; i < next.length; i++) {
			if (next[i] < max[i]) {
				next[i]++;
				return;
			} else {
				next[i] = min[i];
			}
		}

		finished = true;
	}

	@Override
	public Iterator<int[]> iterator() {
		reset();
		return this;
	}

	@Override
	public int[] next() {
		if (inc) inc();
		if (!hasNext()) return null;

		inc = true;
		return next;
	}

	@Override
	public void remove() {

	}

	void reset() {
		inc = false;
		finished = false;
		for (int i = 0; i < next.length; i++)
			next[i] = min[i];
	}

	public void set(int idx, int min, int max) {
		if (max <= min) throw new RuntimeException("Cannot initialize 'max' less or equal than 'min'. min = " + min + ", max = " + max);
		this.min[idx] = min;
		this.max[idx] = max;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[ ");
		for (int i = 0; i < next.length; i++)
			sb.append((i > 0 ? ", " : "") + next[i]);
		sb.append(" ]");

		return sb.toString();
	}

}
