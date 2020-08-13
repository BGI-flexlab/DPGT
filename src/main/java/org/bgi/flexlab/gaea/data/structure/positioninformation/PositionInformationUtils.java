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
package org.bgi.flexlab.gaea.data.structure.positioninformation;

import java.util.ArrayList;

public class PositionInformationUtils<T> {
	protected ArrayList<T> info;
	
	public PositionInformationUtils() {
		info = new ArrayList<T>();
	}
	
	public PositionInformationUtils(int windowSize) {
		info = new ArrayList<T>(windowSize);
	}
	
	public T get(int i) {
		return info.get(i);
	}
	
	public void add(T t) {
		info.add(t);
	}
	

}
