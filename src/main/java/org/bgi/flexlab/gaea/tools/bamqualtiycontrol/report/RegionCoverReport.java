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
package org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report;

import java.util.Arrays;

public class RegionCoverReport {
	private int[] depth;
	
	public RegionCoverReport() {
		depth = new int[1000 + 1];
		init();
	}
	
	public RegionCoverReport(int windowSize) {
		depth = new int[windowSize + 1];
		init();
	}
	
	private void init() {
		Arrays.fill(depth, 0);
	}
	
	public void add(int depth) {
		if(depth < this.depth.length) {
			this.depth[depth]++;
		}
	}
	
	public String toReducerString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Region Depth:\n");
		sb.append(depth[0]);
		for(int i = 1; i < depth.length; i++) {
			sb.append("\t");
			sb.append(depth[i]);
		}
		sb.append("\n");
		return sb.toString();
	}
	
	//name = "RMDUP Rgion Depth"
	public String toReducerString(String name) {
		StringBuilder sb = new StringBuilder();
		sb.append(name + ":\n");
		sb.append(depth[0]);
		for(int i = 1; i < depth.length; i++) {
			sb.append("\t");
			sb.append(depth[i]);
		}
		sb.append("\n");
		return sb.toString();
	}
	
	public void parseReducerOutput(String line) {
		String[] lineSplits = line.split("\t");
		for(int i = 0; i < lineSplits.length; i++) {
			depth[i] += Integer.parseInt(lineSplits[i]);
		}
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < depth.length; i++) {
			sb.append(i);
			sb.append("\t");
			sb.append(depth[i]);
			sb.append("\n");
		}
		return sb.toString();
	}
	
	public int[] getDepthArray() {
	    return depth;
	}
}
