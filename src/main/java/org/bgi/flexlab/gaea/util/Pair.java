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
 *******************************************************************************/
package org.bgi.flexlab.gaea.util;

public class Pair<X, Y> {
	// declare public, STL-style for easier and more efficient access:
	public X first;
	public Y second;

	public Pair(X x, Y y) {
		first = x;
		second = y;
	}

	public void set(X x, Y y) {
		first = x;
		second = y;
	}

	/**
	 * Utility method to make creating pairs easier by skipping typing out the types
	 * @param one the first item
	 * @param two the second item
	 * @param <T1> type of the first item
	 * @param <T2> type of the second item
	 * @return a pair containing both items
	 */
	public static<T1, T2> Pair<T1, T2> create(T1 one, T2 two) {
		return new Pair<>(one, two);
	}

	public X getFirst() {
		return first;
	}

	public Y getSecond() {
		return second;
	}

	@SuppressWarnings("rawtypes")
	@Override
	public boolean equals(Object object) {
		if (object == null){
			return false;
		}
		if (!(object instanceof Pair)){
			return false;
		}

		Pair other = (Pair) object;

		// Check to see whether one is null but not the other.
		if (this.first == null && other.first != null){
			return false;
		}
		if (this.second == null && other.second != null){
			return false;
		}

		// Check to see whether the values are equal.
		if (this.first != null && !this.first.equals(other.first)){
			return false;
		}
		if (this.second != null && !this.second.equals(other.second)){
			return false;
		}

		return true;
	}

	@Override
	public int hashCode() {
		if (second == null && first == null){
			return 0;
		}
		if (second == null){
			return first.hashCode();
		}
		if (first == null){
			return second.hashCode();
		}
		return first.hashCode() ^ second.hashCode();
	}

	public String toString() {
		return first + "," + second;
	}
}