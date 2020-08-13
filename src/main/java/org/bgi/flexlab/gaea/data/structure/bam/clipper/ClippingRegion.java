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
package org.bgi.flexlab.gaea.data.structure.bam.clipper;

public class ClippingRegion {
	/**
	 * contains start position
	 */
	private final int start;
	
	/**
	 * not contains stop position
	 */
	private final int stop;

	public ClippingRegion(int start, int stop) {
		this.start = start;
		this.stop = stop;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public int getLength() {
		return stop - start;
	}
}
