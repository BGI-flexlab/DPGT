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

import org.bgi.flexlab.gaea.data.structure.positioninformation.IntPositionInformation;
import org.bgi.flexlab.gaea.data.structure.positioninformation.depth.PositionDepth;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegion;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegion.Regiondata;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegionStatistic;
import org.bgi.flexlab.gaea.util.ChromosomeUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class SingleRegionReport {
	protected SingleRegion singleReigon;
	protected StringBuffer outputString;
	private double allRegionCoverage = 0;
	private double allRegionFourxCoverage = 0;
	private double allRegionTenxCoverage = 0;
	private double allRegionThirtyxCoverage = 0;
	private double allRegionAverageDepth = 0;
	private double allRegionAverageRmdupDepth = 0;
	private double regionSizeTotal = 0;
	private boolean statPosDepth = false;
	protected Map<SingleRegion.Regiondata, SingleRegionStatistic> result = null;
	
	public SingleRegionReport(SingleRegion singleReigon) {
		this.singleReigon = singleReigon;
		outputString = new StringBuffer();
		result = new ConcurrentHashMap<>();
	}
	
	public String toReducerString() {
		return outputString.toString();
	}

	public String getStatisticString(String chrName, int winStart, int windowSize, PositionDepth dp, String title) {
		int start = 0, end, index = -1;
		end = windowSize - 1;

		IntPositionInformation rmdupDeep = dp.getRMDupPosDepth();
		IntPositionInformation deep = null;
		if(title.equalsIgnoreCase("gender"))
			deep = dp.getGenderPosDepth();
		else
			deep = dp.getNormalPosDepth();

		int i = start;
		while(i <= end) {
			if((index = singleReigon.posInRegion(chrName, i + winStart)) >= 0) {
				break;
			}
			else {
				i++;
			}
		}

		if(index >= 0) {
			while(withinChrAndBin(index, winStart, windowSize, chrName)) {
				int regionStart = singleReigon.getRegion(index).getStart();
				int regionEnd = singleReigon.getRegion(index).getEnd();
				//System.err.println("index:" + index + "\tstart-end:" +regionStart+"-"+regionEnd + "\tposition:"+i +"\t" + winStart);
				if(regionStart >= winStart && regionEnd <= winStart + windowSize -1) {
					if(isStatPosDepth()) {
						outputString.append(title + " part single Region Statistic:\n" + index + ":");
						getPartRegionInfo(rmdupDeep, deep, regionStart - winStart, regionEnd - winStart);
					}else {
						outputString.append(title + " single Region Statistic:\n" + index + ":");
						getWholeRegionInfo(rmdupDeep, deep,regionStart - winStart, regionEnd - winStart);
					}
				} else if(regionStart >= winStart && regionEnd > winStart + windowSize - 1) {
					outputString.append(title + " part single Region Statistic:\n" + index + ":");
					getPartRegionInfo(rmdupDeep, deep, regionStart - winStart, windowSize - 1);
				} else {
					outputString.append(title + " part single Region Statistic:\n" + index + ":");
					if(regionEnd- winStart >= windowSize) {
						getPartRegionInfo(rmdupDeep, deep, 0, windowSize - 1);
					} else {
						getPartRegionInfo(rmdupDeep, deep,0, regionEnd- winStart);
					}
				}
				index++;
			}
		}
		return outputString.toString();
	}
	
	private boolean withinChrAndBin(int index, int winStart, int windowSize, String chrName) {
		return index < singleReigon.getChrInterval(chrName)[1] && 
				singleReigon.getRegion(index).getStart() < winStart + windowSize && 
				singleReigon.getRegion(index).getEnd() >= winStart;
	}
	public Map<Regiondata, SingleRegionStatistic> getResult() {
		return result;
	}
	
	public SingleRegion getRegion() {
		return singleReigon;
	}

	public String getWholeRegionInfo(IntPositionInformation rmdupDeep, IntPositionInformation deep, int start, int end) {
		ArrayList<Integer> deepth = new ArrayList<Integer>();
		ArrayList<Integer> rmdupDeepth = new ArrayList<Integer>();
		SingleRegionStatistic statistic = new SingleRegionStatistic();

		if(start > end) {
			throw new RuntimeException("start:" + start + "is bigger than end:" + end);
		}
//		int n = 0;
		for(int i = start; i <= end; i++) {
			int deepNum = deep.get(i);
			int rmdupDeepNum = rmdupDeep.get(i);
			if(rmdupDeepNum > 0) {
//				n++;
				statistic.addCoverage(rmdupDeepNum);
				statistic.addDepth(deepNum);
				statistic.addRmdupDepth(rmdupDeepNum);
			}
			deepth.add(deepNum);
			rmdupDeepth.add(rmdupDeepNum);
		}
		Collections.sort(deepth);
		Collections.sort(rmdupDeepth);
//		System.out.println("deep Num > 0:" + n);
		outputString.append(statistic.toString());
		outputString.append("\t");
		outputString.append(statistic.calMiddleVaule(deepth));
		outputString.append("\t");
		outputString.append(statistic.calMiddleVaule(rmdupDeepth));
		outputString.append("\n");
		return outputString.toString();
	}

	public String getPartRegionInfo(IntPositionInformation rmdupDeep, IntPositionInformation deep, int start, int end) {
		SingleRegionStatistic statistic = new SingleRegionStatistic();

		if(rmdupDeep.get(start) != 0){
			statistic.addCoverage(rmdupDeep.get(start));
			statistic.addDepth(deep.get(start));
			statistic.addRmdupDepth(rmdupDeep.get(start));
		}
		outputString.append(deep.get(start)).append(",").append(rmdupDeep.get(start));
		for(int i = start + 1; i <= end; i++) {
			outputString.append("\t");
			outputString.append(deep.get(i)).append(",").append(rmdupDeep.get(i));
			if(rmdupDeep.get(i) != 0){
				statistic.addCoverage(rmdupDeep.get(i));
				statistic.addDepth(deep.get(i));
				statistic.addRmdupDepth(rmdupDeep.get(i));
			}
		}
		outputString.append(":");
		outputString.append(statistic.toString());
		outputString.append("\n");

		return outputString.toString();
	}

	public void parseReducerOutput(String line, boolean isPart) {
		String[] splits = line.split(":");
		if(isPart)
			updatePartResult(splits);
		else
			updateResult(splits);
	}

	private void updatePartResult(String[] splits) {
		ArrayList<String> depth = new ArrayList<String>();
		int regionIndex = Integer.parseInt(splits[0]);
		Collections.addAll(depth, splits[1].split("\t"));
		String[] coverSplits = splits[2].split("\t");
		int coverNum = Integer.parseInt(coverSplits[0]);
		int fourxCoverNum = Integer.parseInt(coverSplits[3]);
		int tenxCoverNum = Integer.parseInt(coverSplits[4]);
		int thirtyxCoverNum = Integer.parseInt(coverSplits[5]);
		int depthNum = Integer.parseInt(coverSplits[1]);
		int rmdupDepthNum = Integer.parseInt(coverSplits[2]);

		SingleRegionStatistic statistic;
		if(!result.containsKey(singleReigon.getRegion(regionIndex))) {
			statistic = new SingleRegionStatistic(true);
			statistic.addPart(coverNum, fourxCoverNum, tenxCoverNum,thirtyxCoverNum, depthNum, rmdupDepthNum, depth);
			result.put(singleReigon.getRegion(regionIndex), statistic);
		} else {
			statistic = result.get(singleReigon.getRegion(regionIndex));
			statistic.addPart(coverNum, fourxCoverNum, tenxCoverNum,thirtyxCoverNum, depthNum, rmdupDepthNum, depth);
		}
	}

	private void updateResult(String[] splits) {
		int regionIndex = Integer.parseInt(splits[0]);
		String[] coverSplits = splits[1].split("\t");
		int coverNum = Integer.parseInt(coverSplits[0]);
		int depthNum =  Integer.parseInt(coverSplits[1]);
		int rmdupDepthNum =  Integer.parseInt(coverSplits[2]);
		int fourxCoverNum = Integer.parseInt(coverSplits[3]);
		int tenxCoverNum = Integer.parseInt(coverSplits[4]);
		int thirtyxCoverNum = Integer.parseInt(coverSplits[5]);
		double middleDepth = Double.parseDouble(coverSplits[6]);
		double middleRmdupDepth = Double.parseDouble(coverSplits[7]);

		SingleRegionStatistic statistic = new SingleRegionStatistic(false);
		statistic.add(coverNum, fourxCoverNum, tenxCoverNum,thirtyxCoverNum, depthNum, rmdupDepthNum, middleDepth, middleRmdupDepth);
		result.put(singleReigon.getRegion(regionIndex), statistic);
	}

	public void updateAllRegionAverageDeepth() {
		long depthAll = 0;
		long rmdupDepthAll = 0;
		long coverBaseAll = 0;
		long fourxCoverBaseAll = 0;
		long tenxCoverBaseAll = 0;
		long thirtyxCoverBaseAll = 0;
		regionSizeTotal = 0.0;
		updateResult();
		System.out.println(result.keySet().size());
		for(Regiondata regionData : result.keySet()) {
			SingleRegionStatistic statistic = (SingleRegionStatistic) result.get(regionData);
			String formatChrName = ChromosomeUtils.formatChrName(regionData.getChrName());
			if(!formatChrName.equals("chrx") && !formatChrName.equals("chry") && !formatChrName.equals("chrm")) {
				depthAll += statistic.getDepthNum();
				rmdupDepthAll += statistic.getRmdupDepthNum();
				coverBaseAll += statistic.getCoverBaseNum();
				fourxCoverBaseAll += statistic.getFourxCoverBaseNum();
				tenxCoverBaseAll += statistic.getTenxCoverBaseNum();
				thirtyxCoverBaseAll += statistic.getThirtyxCoverBaseNum();
				regionSizeTotal += regionData.size();
			}
		}
		System.err.println("depthAll:" + depthAll + "\nregionSizeTotal:" + regionSizeTotal);
		allRegionAverageDepth = depthAll / regionSizeTotal;
		allRegionAverageRmdupDepth = rmdupDepthAll / regionSizeTotal;
		allRegionCoverage = coverBaseAll / regionSizeTotal;
		allRegionFourxCoverage = fourxCoverBaseAll / regionSizeTotal;
		allRegionTenxCoverage = tenxCoverBaseAll / regionSizeTotal;
		allRegionThirtyxCoverage = thirtyxCoverBaseAll / regionSizeTotal;
		System.err.println("allRegionAverageDeepth:" + allRegionAverageDepth);
	}

	private void updateResult() {
		for(Regiondata regionData : singleReigon.getRegions()) {
			if(!result.containsKey(regionData)) {
				//System.err.println("no region:" + regionData.getNameString());
				SingleRegionStatistic statistic = new SingleRegionStatistic(false);
				result.put(regionData, statistic);
			}
		}
	}

	public double getAllRegionAverageDepth() {
		return allRegionAverageDepth;
	}

	public double getAllRegionAverageRmdupDepth() {
		return allRegionAverageRmdupDepth;
	}

	public double getRegionSizeTotal() {
		return regionSizeTotal;
	}

	public double getAllRegionCoverage() {
		return allRegionCoverage;
	}

	public double getAllRegionFourxCoverage() {
		return allRegionFourxCoverage;
	}

	public double getAllRegionTenxCoverage() {
		return allRegionTenxCoverage;
	}

	public double getAllRegionThirtyxCoverage() {
		return allRegionThirtyxCoverage;
	}

	public boolean isStatPosDepth() {
		return statPosDepth;
	}

	public void setStatPosDepth(boolean statPosDepth) {
		this.statPosDepth = statPosDepth;
	}
}
