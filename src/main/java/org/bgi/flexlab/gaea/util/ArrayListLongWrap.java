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

import java.util.ArrayList;

public class ArrayListLongWrap {
	private ArrayList<Long> data;

	public ArrayListLongWrap() {
		data = new ArrayList<Long>();
	}

	// index start from 0
	public void add(int index, long value) {
		if (index >= data.size()) {
			for (int i = data.size(); i <= index - 1; i++) {
				data.add((long) 0);
			}
			data.add(value);
		} else {
			data.set(index, data.get(index) + value);
		}
	}

	public int size() {
		return data.size();
	}

	public String toString() {
		if (data.size() == 0)
			return "";
		StringBuilder sb = new StringBuilder();
		sb.append(data.get(0));
		for (int i = 1; i < data.size(); i++) {
			sb.append("\t");
			sb.append(data.get(i));
		}

		return sb.toString();
	}

	public void readString(String line) {
		if (line.trim().length() == 0)
			return;
		String[] LineSplits = line.split("\t");
		int i = 0;
		for (String split : LineSplits) {
			add(i, Long.parseLong(split));
			i++;
		}
	}

	public ArrayList<Long> get() {
		return this.data;
	}

	public void add(ArrayListLongWrap other) {
		ArrayList<Long> _other = other.get();
		int size = _other.size();

		for (int i = 0; i < size; i++)
			add(i, _other.get(i));
	}
}
