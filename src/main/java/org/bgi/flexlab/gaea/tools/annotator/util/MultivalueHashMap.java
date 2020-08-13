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

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * A Hash that can hold multiple values for each key
 *
 * @author pcingola
 */
public class MultivalueHashMap<K, V> extends HashMap<K, LinkedList<V>> {

	private static final long serialVersionUID = -8860279543165990227L;

	public MultivalueHashMap() {
		super();
	}

	/**
	 * Add multiple values
	 */
	public void add(K key, Collection<V> values) {
		getOrCreate(key).addAll(values); // Add all to the list
	}

	/**
	 * Add a single value
	 */
	public void add(K key, V value) {
		getOrCreate(key).add(value); // Add to the list
	}

	/**
	 * Get a list of values (or create it if not available)
	 */
	public List<V> getOrCreate(K key) {
		// Get list
		LinkedList<V> list = get(key);
		if (list == null) { // No list? Create one
			list = new LinkedList<V>();
			put(key, list);
		}
		return list;
	}
}
