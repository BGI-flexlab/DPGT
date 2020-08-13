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
package org.bgi.flexlab.gaea.tools.annotator.db;

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class Results extends HashMap<String, LinkedList<HashMap<String,String>>> {

	private static final long serialVersionUID = -8860279543165990227L;

	private Map<String, HashMap<String, String>> mergeResults = new HashMap<>();

	public Results() {
		super();
	}

	/**
	 * Add multiple values
	 */
	public void add(String key, Collection<HashMap<String,String>> values) {
		getOrCreate(key).addAll(values); // Add all to the list
	}

	/**
	 * Add a single value
	 */
	public void add(String key, HashMap<String,String> value) {
		getOrCreate(key).add(value); // Add to the list
	}

	/**
	 * Get a list of values (or create it if not available)
	 */
	public List<HashMap<String,String>> getOrCreate(String key) {
		// Get list
		LinkedList<HashMap<String,String>> list = get(key);
		if (list == null) { // No list? Create one
			list = new LinkedList<HashMap<String,String>>();
			put(key, list);
		}
		return list;
	}

	/**
	 * Get a mreged result, key is: alt or gene
	 */
	public HashMap<String,String> getMergeResult(String key) {
		// Get list
		LinkedList<HashMap<String,String>> list = get(key);
		if (list == null)
			return null;
		if(mergeResults.containsKey(key))
			return mergeResults.get(key);
		else {
			HashMap<String,String> r = mergeResult(list);
			mergeResults.put(key, r);
			return r;
		}
	}

	private HashMap<String, String> mergeResult(
			LinkedList<HashMap<String, String>> resultList) {
		if (resultList == null || resultList.isEmpty()) return null;

		HashMap<String, String> result = resultList.poll();
		HashMap<String, String> temp;
		while ((temp = resultList.poll()) != null){
			for (Entry<String, String> entry : temp.entrySet()) {
				String key = entry.getKey();
				String value = entry.getValue();
				if(result.putIfAbsent(key, value) != null) {
					if (!result.get(key).equals(value)) {
						String mergeValue = result.get(key) + "&" + value;
						result.put(key, mergeValue);
					}
				}
			}
		}
		return result;
	}
}
