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
import java.util.Map;

/**
 * Tuple: A pair of objects
 * @author pablocingolani
 *
 * @param <A>
 * @param <B>
 */
@SuppressWarnings("rawtypes")
public class Tuple<A, B> {
	public final A first;
	public final B second;

	public Tuple(A first, B second) {
		super();
		this.first = first;
		this.second = second;
	}

	@Override
	public boolean equals(Object other) {
		if (other instanceof Tuple) {
			Tuple otherPair = (Tuple) other;
			return ((this.first == otherPair.first || (this.first != null && otherPair.first != null && this.first.equals(otherPair.first))) && (this.second == otherPair.second || (this.second != null && otherPair.second != null && this.second.equals(otherPair.second))));
		}
		return false;
	}

	public A getFirst() {
		return first;
	}

	public B getSecond() {
		return second;
	}

	@Override
	public int hashCode() {
		int hashFirst = first != null ? first.hashCode() : 0;
		int hashSecond = second != null ? second.hashCode() : 0;

		return (hashFirst + hashSecond) * hashSecond + hashFirst;
	}

	public static Tuple<String, String> transMapToTuple(Map map){
		return transMapToTuple(map, ";");
	}

	public static Tuple<String, String> transMapToTuple(Map map, String sep){
		java.util.Map.Entry entry;
		StringBuilder k = new StringBuilder();
		StringBuilder v = new StringBuilder();
		for(Iterator iterator = map.entrySet().iterator(); iterator.hasNext();)
		{
			entry = (java.util.Map.Entry)iterator.next();
			k.append(entry.getKey().toString()).append (iterator.hasNext() ? ";" : "");
			v.append(entry.getValue().toString()).append (iterator.hasNext() ? ";" : "");
		}
		return new Tuple<>(k.toString(), v.toString());
	}

	@Override
	public String toString() {
		return "(" + first + ", " + second + ")";
	}
}