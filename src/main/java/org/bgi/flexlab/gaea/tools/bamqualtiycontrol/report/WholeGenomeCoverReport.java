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

import org.bgi.flexlab.gaea.data.structure.positioninformation.depth.PositionDepth;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.BaseCounter;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.BaseType;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.Depth;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.DepthType;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.Interval;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.Tracker;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.Tracker.BaseTracker;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;



public class WholeGenomeCoverReport{

	private BaseTracker bTracker;

	private ChromosomeInformationShare chrInfo;
		
	private static Map<String, WholeGenomeCoverReport> coverReports = new ConcurrentHashMap<>();
	
	public WholeGenomeCoverReport(ChromosomeInformationShare chrInfo) {
		bTracker = new Tracker.BaseTracker();
		register();
		this.chrInfo = chrInfo;
	}
	
	public void constructDepthReport(PositionDepth deep, int i, int pos) {
		int depth = deep.getPosDepth(i);
		int rmdupDepth = deep.getRMDupPosDepth(i);
		if(depth != 0) {
			bTracker.setTrackerAttribute(Interval.WHOLEGENOME, DepthType.NORMAL.setDepth(depth), DepthType.WITHOUT_PCR.setDepth(rmdupDepth));
			bTracker.setTrackerAttribute(BaseType.COVERED);
			if(rmdupDepth >= 4){
				bTracker.setTrackerAttribute(BaseType.FOURXCOVERED);
				if(chrInfo.getBase(pos) != 'N')
					bTracker.setTrackerAttribute(BaseType.FOURXNONNCOVERED);
			}
			if(rmdupDepth >= 10){
				bTracker.setTrackerAttribute(BaseType.TENXCOVERED);
				if(chrInfo.getBase(pos) != 'N')
					bTracker.setTrackerAttribute(BaseType.TENXNONNCOVERED);
			}
			if(rmdupDepth >= 30){
				bTracker.setTrackerAttribute(BaseType.THIRTYXCOVERED);
				if(chrInfo.getBase(pos) != 'N')
					bTracker.setTrackerAttribute(BaseType.THIRTYXNONNCOVERED);
			}
			if(chrInfo.getBase(pos) != 'N')
				bTracker.setTrackerAttribute(BaseType.NONNCOVERED);
			if(deep.hasIndelReads(i))
				bTracker.setTrackerAttribute(BaseType.INDELREF);
			
			if(deep.hasMismatchReads(i)) 
				bTracker.setTrackerAttribute(BaseType.MISMATCHREF);				
		} else {
			if(deep.isDeletionBaseWithNoConver(i)) {
				bTracker.setTrackerAttribute(BaseType.COVERED);
//				bTracker.setTrackerAttribute(BaseType.FOURX);
//				TODO  推测delelion附近深度
				if(chrInfo.getBase(pos) != 'N')
					bTracker.setTrackerAttribute(BaseType.NONNCOVERED);
			}
		}
	}
	
	public String toReducerString(String chrName) {
		StringBuffer coverString = new StringBuffer();
		coverString.append("Cover Information:\n");
		coverString.append(chrName);
		coverString.append("\t");
		for(Entry<String, BaseCounter> counter : bTracker.getCounterMap().entrySet()) {
			coverString.append(counter.getKey());
			coverString.append(" ");
			coverString.append(counter.getValue().getProperty());
			coverString.append("\t");
		}
		coverString.append("\n");
		
		return coverString.toString();
	}
	
	public String toString() {
		DecimalFormat df = new DecimalFormat("0.000");
		df.setRoundingMode(RoundingMode.HALF_UP);
		
		StringBuffer coverString = new StringBuffer();
		coverString.append("chromsome:\t");
		coverString.append(chrInfo.getChromosomeName());
		coverString.append("\nCoverage (>=1x):\t");
		coverString.append(df.format(getCoverage()));
		coverString.append("%\nCoverage (>=4x):\t");
		coverString.append(df.format(getFourxCoverage()));
		coverString.append("%\nCoverage (>=10x):\t");
		coverString.append(df.format(getTenxCoverage()));
		coverString.append("%\nCoverage (>=30x):\t");
		coverString.append(df.format(getThirtyxCoverage()));
		coverString.append("%\nMean Depth:\t");
		coverString.append(df.format(getMeanDepth()));
		coverString.append("\nMean Rmdup Depth:\t");
		coverString.append(df.format(getMeanRmdupDepth()));
		coverString.append("\nNonN Coverage (>=1x):\t");
		coverString.append(df.format(getNonNbaseCoverage()));
		coverString.append("%\nNonN Coverage (>=4x):\t");
		coverString.append(df.format(getFourxNonNCoverage()));
		coverString.append("%\nNonN Coverage (>=10x):\t");
		coverString.append(df.format(getTenxNonNCoverage()));
		coverString.append("%\nNonN Coverage (>=30x):\t");
		coverString.append(df.format(getThirtyxNonNCoverage()));
		coverString.append("%\nNonN Mean Depth:\t");
		coverString.append(df.format(getNonNMeanDepth()));
		coverString.append("\nNonN Mean Rmdup Depth:\t");
		coverString.append(df.format(getNonNMeanRmdupDepth()));
		coverString.append("\nrate of position according to reference that have at least one indel reads support:\t");
		coverString.append(df.format(getRateOf(BaseType.INDELREF)));
		coverString.append("%\nrate of position according to reference that have at least one mismatch reads support:\t");
		coverString.append(df.format(getRateOf(BaseType.MISMATCHREF)));
		coverString.append("%\n\n");
		
		return coverString.toString();
	}
	
	public void parse(String keyValue, ReferenceShare genome) {
		String key = keyValue.split(" ")[0];
		String value = keyValue.split(" ")[1];
		BaseCounter bCounter = null;
		if((bCounter = getBasetTracker().getCounterMap().get(key)) != null) {
			if(!key.contains(Depth.TOTALDEPTH.toString())) {
				bCounter.setBaseCount(Long.parseLong(value));
			} else if(key.contains(DepthType.NORMAL.toString())) {
				bCounter.setTotalDepth(Long.parseLong(value));
			}else {
				bCounter.setTotalDepthWithoutPCRDup(Long.parseLong(value));
			}
		} else {
			throw new RuntimeException("Can not idenity counter with name " + key);
		}
	}
	
	public BaseTracker getBasetTracker() {
		return bTracker;
	}

	public double getDepth() {
		return bTracker.getProperty(Interval.WHOLEGENOME, Depth.TOTALDEPTH, DepthType.NORMAL);
	}

	public double getRmdupDepth() {
		return bTracker.getProperty(Interval.WHOLEGENOME, Depth.TOTALDEPTH, DepthType.WITHOUT_PCR);
	}

	public long getCoverBaseNum() {
		return bTracker.getProperty(BaseType.COVERED);
	}

	public long getFourXCoverBaseNum() {
		return bTracker.getProperty(BaseType.FOURXCOVERED);
	}

	public long getTenXCoverBaseNum() {
		return bTracker.getProperty(BaseType.TENXCOVERED);
	}

	public long getTenXNonNCoverBaseNum() {
		return bTracker.getProperty(BaseType.TENXNONNCOVERED);
	}

	public long getThirtyXCoverBaseNum() {
		return bTracker.getProperty(BaseType.THIRTYXCOVERED);
	}

	public long getThirtyXNonnCoverBaseNum() {
		return bTracker.getProperty(BaseType.THIRTYXNONNCOVERED);
	}

	public long getNonNCoverBaseNum() {
		return bTracker.getProperty(BaseType.NONNCOVERED);
	}

	public long getFourXNonNCoverBaseNum() {
		return bTracker.getProperty(BaseType.FOURXNONNCOVERED);
	}


	public double getCoverage() {
		return (100 * (getCoverBaseNum()/(double)chrInfo.getLength()));
	}

	public double getFourxCoverage() {
		return (100 * (getFourXCoverBaseNum()/(double)chrInfo.getLength()));
	}

	public double getFourxNonNCoverage() {
		return (100 * (getFourXNonNCoverBaseNum()/(double)chrInfo.getNonNbaselength()));
	}

	public double getTenxNonNCoverage() {
		return (100 * (getTenXNonNCoverBaseNum()/(double)chrInfo.getNonNbaselength()));
	}

	public double getThirtyxNonNCoverage() {
		return (100 * (getThirtyXNonnCoverBaseNum()/(double)chrInfo.getNonNbaselength()));
	}

	public double getTenxCoverage() {
		return (100 * (getTenXCoverBaseNum()/(double)chrInfo.getLength()));
	}

	public double getThirtyxCoverage() {
		return (100 * (getThirtyXCoverBaseNum()/(double)chrInfo.getLength()));
	}

	public double getNonNbaseCoverage() {
		return (100 * (getNonNCoverBaseNum()/(double)chrInfo.getNonNbaselength()));
	}

	public double getNonNMeanDepth() {
		return (getDepth()/(double)bTracker.getProperty(BaseType.NONNCOVERED));
	}

	public double getNonNMeanRmdupDepth() {
		return (getRmdupDepth()/(double)bTracker.getProperty(BaseType.NONNCOVERED));
	}
	
	public double getMeanDepth() {
		return (getDepth()/(double)bTracker.getProperty(BaseType.COVERED));
	}

	public double getMeanRmdupDepth() {
		return (getRmdupDepth()/(double)bTracker.getProperty(BaseType.COVERED));
	}

	public double getRateOf(BaseType type) {
		return (100 * (bTracker.getProperty(type)/(double)chrInfo.getLength()));
	}
	
	public void register() {
		bTracker.register(createBaseCounters());
	}

	public List<BaseCounter> createBaseCounters() {
		List<BaseCounter> counters = new ArrayList<>();
		Collections.addAll(counters, new BaseCounter(Interval.WHOLEGENOME, Depth.TOTALDEPTH, DepthType.NORMAL),
				                     new BaseCounter(Interval.WHOLEGENOME, Depth.TOTALDEPTH, DepthType.WITHOUT_PCR),
									 new BaseCounter(BaseType.FOURXCOVERED),
									 new BaseCounter(BaseType.FOURXNONNCOVERED),
									 new BaseCounter(BaseType.TENXCOVERED),
									 new BaseCounter(BaseType.TENXNONNCOVERED),
									 new BaseCounter(BaseType.THIRTYXCOVERED),
									 new BaseCounter(BaseType.THIRTYXNONNCOVERED),
									 new BaseCounter(BaseType.NONNCOVERED),
									 new BaseCounter(BaseType.COVERED),
									 new BaseCounter(BaseType.INDELREF),
									 new BaseCounter(BaseType.MISMATCHREF));

		return counters;
	}
	
	public static Map<String, WholeGenomeCoverReport> getCoverReports() {
		return coverReports;
	}
	
	public static WholeGenomeCoverReport getCoverReport(ChromosomeInformationShare chrInfo) {
		if(!coverReports.containsKey(chrInfo.getChromosomeName())) {
			WholeGenomeCoverReport coverInfo = new WholeGenomeCoverReport(chrInfo);
			coverReports.put(chrInfo.getChromosomeName(), coverInfo);
		}
		return coverReports.get(chrInfo.getChromosomeName());
	}
	
	public static WholeGenomeCoverReport getCoverReport(String chrName) {
		return coverReports.get(chrName);
	}
	
	public static void addCoverReport(ChromosomeInformationShare chrInfo) {
		WholeGenomeCoverReport coverInfo = new WholeGenomeCoverReport(chrInfo);
		coverReports.put(chrInfo.getChromosomeName(), coverInfo);
	}

	public ChromosomeInformationShare getChrInfo() {
		return chrInfo;
	}

}
