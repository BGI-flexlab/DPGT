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
package org.bgi.flexlab.gaea.data.structure.region;

import org.bgi.flexlab.gaea.data.structure.region.SingleRegion.Regiondata;
import org.bgi.flexlab.gaea.util.MathUtils;

public class RegionChromosomeInformation {
	private int size;
	private int coverageNum;
	private long depthNum;
	
	public RegionChromosomeInformation() {
		size = 0;
		coverageNum = 0;
		depthNum = 0;
	}
	
	public void add(Regiondata regiondata, SingleRegionStatistic singleRegionInfo) {
		size += regiondata.size();
		coverageNum += singleRegionInfo.getCoverBaseNum();
		depthNum += singleRegionInfo.getDepthNum();
	}

	public static String toTitleString() {
		return "Chr\tRegionSize\tCoverage\tAverageDepth\n";
	}

	public String toString(String chrName) {
		StringBuilder sb = new StringBuilder();
		sb.append(chrName);
		sb.append("\t");
		sb.append(size);
		sb.append("\t");
		sb.append(MathUtils.doubleformat.format(getCoverage()));
		sb.append("\t");
		sb.append(MathUtils.doubleformat.format(getAverageDepth()));
		sb.append("\n");
		
		return sb.toString();
	}

	/**
	 * @return the averageDepth
	 */
	public double getAverageDepth() {
		return depthNum / (double) coverageNum;
	}

	/**
	 * @return the coverage
	 */
	public double getCoverage() {
		return coverageNum / (double) size;
	}


	/**
	 * @return the size
	 */
	public int getSize() {
		return size;
	}
}
