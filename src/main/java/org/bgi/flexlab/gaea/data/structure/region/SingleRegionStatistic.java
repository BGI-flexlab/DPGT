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

import org.bgi.flexlab.gaea.data.exception.BAMQCException;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.util.StatsUtils;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;

public class SingleRegionStatistic {
	private double refGCrate;
	private int coverBaseNum = 0;
	private int fourxCoverBaseNum = 0;
	private int tenxCoverBaseNum = 0;
	private int thirtyxCoverBaseNum = 0;
	private int depthNum = 0;
	private int rmdupDepthNum = 0;
	private double middleDepth = 0;
	private double middleRmdupDepth = 0;
	private ArrayList<Integer> depth = null;
	private ArrayList<Integer> rmdupDepth = null;

	public SingleRegionStatistic() {}

	public SingleRegionStatistic(boolean isPart) {
		if(isPart && depth == null) {
			depth = new ArrayList<Integer>();
			rmdupDepth = new ArrayList<Integer>();
		}
	}

	public SingleRegionStatistic(int depthNum, int coveragePosNum) {
		this.depthNum = depthNum;
		this.coverBaseNum = coveragePosNum;
		refGCrate = 0;
	}

	public double calMiddleVaule(ArrayList<Integer> deepth) {
		double middleValue;
		if(deepth.size() == 0) {
			middleValue = 0;
			return middleValue;
		}
		if(deepth.size() % 2 == 0) {
			middleValue = (deepth.get(deepth.size() >> 1) + deepth.get((deepth.size() >> 1) - 1)) / (double)2;
		} else {
			middleValue = deepth.get(deepth.size() >> 1);
		}
		return middleValue;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(coverBaseNum);
		sb.append("\t");
		sb.append(depthNum);
		sb.append("\t");
		sb.append(rmdupDepthNum);
		sb.append("\t");
		sb.append(fourxCoverBaseNum);
		sb.append("\t");
		sb.append(tenxCoverBaseNum);
		sb.append("\t");
		sb.append(thirtyxCoverBaseNum);

		return sb.toString();
	}
	
	public void addCoverage(int depthNum) {
		coverBaseNum++;
		if(depthNum >= 4)
			fourxCoverBaseNum++;
		if(depthNum >= 10)
			tenxCoverBaseNum++;
		if(depthNum >= 30)
			thirtyxCoverBaseNum++;
	}
	
	public void addDepth(int depth) {
		depthNum += depth;
	}

	public void addRmdupDepth(int rmdupDepth) {
		rmdupDepthNum += rmdupDepth;
	}

	public static String toReportTitleString() {
		return "#Chr\tStart\tEnd\tAve depth\tMedian\tAve depth(rmdup)\tMedian(rmdup)\tCoverage(>=1x)%\tCoverage(>=4x)%\tCoverage(>=10x)%\tCoverage(>=30x)%\tnormalizedDepth";
	}

	public String toReportString(SingleRegion.Regiondata regiondata, double allRegionAverageDepth) {
		double averageDepth = getDepth(regiondata);
		double averageRmdupDepth = getRmdupDepth(regiondata);
		double normalizedDepth = averageDepth / (double)allRegionAverageDepth;
		if(depth != null) {
			while(depth.size() < regiondata.size()) {
				depth.add(0);
			}
			Collections.sort(depth);
			middleDepth = calMiddleVaule(depth);
		}
		if(rmdupDepth != null) {
			while(rmdupDepth.size() < regiondata.size()) {
				rmdupDepth.add(0);
			}
			Collections.sort(rmdupDepth);
			middleRmdupDepth = calMiddleVaule(rmdupDepth);
		}

		DecimalFormat df = new DecimalFormat("0.0000");
		df.setRoundingMode(RoundingMode.HALF_UP);
		StringBuilder outString = new StringBuilder();
		outString.append(regiondata.getChrName());
		outString.append("\t");
		outString.append(regiondata.getStart());
		outString.append("\t");
		outString.append(regiondata.getEnd());
		outString.append("\t");
		outString.append(df.format(averageDepth));
		outString.append("\t");
		outString.append(df.format(middleDepth));
		outString.append("\t");
		outString.append(df.format(averageRmdupDepth));
		outString.append("\t");
		outString.append(df.format(middleRmdupDepth));
		outString.append("\t");
		outString.append(StatsUtils.perc(coverBaseNum, regiondata.size()));
		outString.append("\t");
		outString.append(StatsUtils.perc(fourxCoverBaseNum, regiondata.size()));
		outString.append("\t");
		outString.append(StatsUtils.perc(tenxCoverBaseNum, regiondata.size()));
		outString.append("\t");
		outString.append(StatsUtils.perc(thirtyxCoverBaseNum, regiondata.size()));
		outString.append("\t");
		outString.append(df.format(normalizedDepth));
		outString.append("\n");

		return outString.toString();
	}

	public void add(String line) {
		String[] lineSplits = line.split("\t");
		if(lineSplits.length != 2) {
			throw new BAMQCException.WrongNumOfColException(2);
		}
		depthNum += Integer.parseInt(lineSplits[0]);
		coverBaseNum += Integer.parseInt(lineSplits[1]);
		/*for(int i = 0; i < 4; i++) {
			ATCG[i] += Short.parseShort(lineSplits[i + 2]);
		}*/
	}

//	public void add(int coverNum, int depthNum, int rmdupDepthNum, double middleDepth, double middleRmdupDepth) {
//		this.coverBaseNum = coverNum;
//		this.depthNum = depthNum;
//		this.rmdupDepthNum = rmdupDepthNum;
//		this.middleDepth = middleDepth;
//		this.middleRmdupDepth = middleRmdupDepth;
//	}
//
//	public void add(int coverNum, int fourxCoverNum, int depthNum, int rmdupDepthNum, double middleDepth, double middleRmdupDepth) {
//		this.coverBaseNum = coverNum;
//		this.fourxCoverBaseNum = fourxCoverNum;
//		this.depthNum = depthNum;
//		this.rmdupDepthNum = rmdupDepthNum;
//		this.middleDepth = middleDepth;
//		this.middleRmdupDepth = middleRmdupDepth;
//	}

	public void add(int coverNum, int fourxCoverNum, int tenxCoverNum, int thirtyxCoverNum,int depthNum, int rmdupDepthNum, double middleDepth, double middleRmdupDepth) {
		this.coverBaseNum = coverNum;
		this.fourxCoverBaseNum = fourxCoverNum;
		this.tenxCoverBaseNum = tenxCoverNum;
		this.thirtyxCoverBaseNum = thirtyxCoverNum;
		this.depthNum = depthNum;
		this.rmdupDepthNum = rmdupDepthNum;
		this.middleDepth = middleDepth;
		this.middleRmdupDepth = middleRmdupDepth;
	}

	public void addPart(int coverNum, int depthNum, int rmdupDepthNum, ArrayList<String> depthPartList) {
		this.coverBaseNum += coverNum;
		this.depthNum += depthNum;
		this.rmdupDepthNum += rmdupDepthNum;

		// depthStr : depth,rmdupDepth
		for(String depthStr : depthPartList) {
			String[] depths = depthStr.split(",");
			depth.add(Integer.valueOf(depths[0]));
			rmdupDepth.add(Integer.valueOf(depths[1]));
		}
	}

	public void addPart(int coverNum, int fourxCoverNum, int tenxCoverNum, int thirtyxCoverNum,int depthNum, int rmdupDepthNum, ArrayList<String> depthPartList) {
		this.coverBaseNum += coverNum;
		this.fourxCoverBaseNum += fourxCoverNum;
		this.tenxCoverBaseNum += tenxCoverNum;
		this.thirtyxCoverBaseNum += thirtyxCoverNum;
		this.depthNum += depthNum;
		this.rmdupDepthNum += rmdupDepthNum;

		// depthStr : depth,rmdupDepth
		for(String depthStr : depthPartList) {
			String[] depths = depthStr.split(",");
			depth.add(Integer.valueOf(depths[0]));
			rmdupDepth.add(Integer.valueOf(depths[1]));
		}
	}

	public void calRefGCrate(ChromosomeInformationShare chrInfo, SingleRegion.Regiondata regionData, int extendSize) {
		int CGbaseNum = 0;
		int ATCGbaseNum = 0;
		if(extendSize == 0) {
			if(regionData.size() < 200) {
				extendSize = (200 - regionData.size()) / 2 + 1;
			}
		}
		for(char base : chrInfo.getGA4GHBaseSequence(regionData.getStart() - extendSize > 0 ? regionData.getStart() - extendSize : 0, regionData.getEnd() + extendSize).toCharArray()) {
			if(base == 'C' || base == 'G' || base == 'c' || base == 'g') {
				CGbaseNum++;
			}
			if(base != 'N' && base != 'n') {
				ATCGbaseNum++;
			}
		}
		refGCrate = CGbaseNum / (double) ATCGbaseNum;
	}

	public double getAverageDepth() {
		if(coverBaseNum == 0) {
			return 0;
		}
		return depthNum / (double) coverBaseNum;
	}

	/**
	 * @return the depthNum
	 */
	public long getDepthNum() {
		return depthNum;
	}

	public long getRmdupDepthNum() {
		return rmdupDepthNum;
	}

	public int getCoverBaseNum() {
		return coverBaseNum;
	}

	public int getFourxCoverBaseNum() {
		return fourxCoverBaseNum;
	}

	public int getTenxCoverBaseNum() {
		return tenxCoverBaseNum;
	}

	public int getThirtyxCoverBaseNum() {
		return thirtyxCoverBaseNum;
	}

	public double getDepth(SingleRegion.Regiondata regiondata) {
		return getDepthNum() /(double) regiondata.size();
	}

	private double getRmdupDepth(SingleRegion.Regiondata regiondata) {
		return getRmdupDepthNum() /(double) regiondata.size();
	}

	public int getPosDepth(int index){
		if(depth == null)
			return 0;
		return depth.get(index);
	}

	public int getPosRmdupDepth(int index){
		if(depth == null)
			return 0;
		return rmdupDepth.get(index);
	}

	public double getRefGCrate() {
		return refGCrate;
	}

	public double getMiddleDepth() {
		return middleDepth;
	}

	public double getMiddleRmdupDepth() {
		return middleRmdupDepth;
	}
}
