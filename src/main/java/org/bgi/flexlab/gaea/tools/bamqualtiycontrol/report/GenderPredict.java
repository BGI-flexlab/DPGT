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

import org.bgi.flexlab.gaea.data.exception.BAMQCException;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegion;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegion.Regiondata;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegionStatistic;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.RegionReport.Sex;
import org.bgi.flexlab.gaea.util.Lowess;

import java.util.*;
import java.util.Map.Entry;


public class GenderPredict {
	private Map<Integer, LowessData> lowessData;
	private List<Double> recaltable;
	ArrayList<Double> GCFinal;
	private Map<String, Integer[]> chrIndex = new HashMap<String, Integer[]>();
	private static double span = 0.001;
	private double totalMeanDepth;
	private int totalRegionNum = 0;
	
	public GenderPredict(Map<Regiondata, SingleRegionStatistic> result, SingleRegion sg, double[] filterVale) {
		lowessData = new HashMap<Integer, GenderPredict.LowessData>();
		GCFinal = new ArrayList<Double>();
		int i = 0;
		String chrName = "";
		Integer[] startEnd = new Integer[2];
		//int totalRegionNum = 0;
		double totalRegionDepth = 0;
		int coverageFilterNum = 0;
				
		//Regiondata rdlast = null;
		for(Regiondata rd : sg.getRegions()) {
			//rdlast = rd;
			SingleRegionStatistic bsrs = result.get(rd);
			if(bsrs == null)
				continue;
			if(bsrs.getAverageDepth() >= filterVale[1] || bsrs.getAverageDepth() <= filterVale[0])
				continue;
			if(bsrs.getCoverBaseNum() / rd.size() < 0.5) {
				coverageFilterNum++;
				continue;
			}
			if(rd.size() < 10)
				continue;
			if(Double.isNaN(bsrs.getRefGCrate()) || bsrs.getRefGCrate() < 0.1 || bsrs.getRefGCrate() > 0.9)
				continue;
			
			lowessData.put(i, new LowessData(bsrs.getAverageDepth(), bsrs.getRefGCrate()));
			totalRegionDepth += bsrs.getAverageDepth();
			totalRegionNum++;
			
			if(!rd.getChrName().equals(chrName)) {
				//chrIndex.put(chrName, i);
				if(chrName == "")
					startEnd[0] = i;
				else {
					startEnd[1] = i - 1;
					Integer[] se = new Integer[2];
					se[0] = startEnd[0];
					se[1] = startEnd[1];
					chrIndex.put(chrName, se);
					startEnd[0] = i;
				}
				chrName = rd.getChrName();
			}
			i++;
		}
		startEnd[1] = i - 1;
		chrIndex.put(chrName, startEnd);
		totalMeanDepth = totalRegionDepth/totalRegionNum;
		System.err.println("origional region number:" + sg.getRegions().size() + "\tcoverage filter region number:" + coverageFilterNum);
		System.err.println("total region number:" + totalRegionNum + "\ttoatal mean depth:" + totalMeanDepth);
	}
	
	public Sex predict() {
		if(totalRegionNum < 10) {
			System.err.println("only less than 10 regions are useable. Abort Gender predict!");
			return Sex.unKown;
		}
		
		Map<String, Double> sexChrsDepth = new HashMap<String, Double>();
		double recalChrDepth = 0;
		double normalChromosomeRegionDepth = 0;
		int normalChromosomeRegionNum = 0;
		
		lowess();
		//System.err.println("lowess end");
		Map<String, Double> chrRd = new HashMap<String, Double>();
		for(String chrName : chrIndex.keySet()) {
			Integer[] se = chrIndex.get(chrName);
			double chrDepth = 0;
			double chrRegionNum = 0;
			List<Double> depthTmp = new ArrayList<Double>();
			for(int i = se[0]; i <= se[1]; i++) {
				double depth = getRecalDepth(lowessData.get(i).getGCrate(), lowessData.get(i).getDepth());
				if(depth < 0)
					throw new RuntimeException("depth < 0 means some wrong in function getRecalDepth()");
				//System.err.println("chrName:" + chrName + "\tindex:" + i + "\trecaldepth:" + depth + "\torgionaldepth:" + lowessData.get(i).getDepth());
				depthTmp.add(depth);
			}
			//System.err.println("finished recal");
			Collections.sort(depthTmp);
			//System.err.println("finished sort");
			int filterLow = (int) (depthTmp.size() * 0.05);
			int filterUp = (int) (depthTmp.size() * 0.95) + 1;
			if(filterUp > depthTmp.size())
				filterUp = depthTmp.size();
			for(int i = filterLow; i < filterUp; i++) {
				double depth = depthTmp.get(i);
				if(!isSexChromosome(chrName)) {
					normalChromosomeRegionDepth += depth;
					normalChromosomeRegionNum++;
				} 
				chrDepth += depth;
				chrRegionNum++;
			}
			//System.err.println("finished filter");
			chrRd.put(chrName, chrDepth/chrRegionNum);
			System.err.println(chrName + "\t" + chrRegionNum + "\t" + chrDepth + "\t" + chrDepth/chrRegionNum);
			if(isSexChromosome(chrName)) {
				if(chrRegionNum > 5)
					sexChrsDepth.put(chrName, chrDepth/chrRegionNum);
				else
					System.err.println("drop chr:" + chrName + " depth due to number of region is below 5.");
			}
		}
		recalChrDepth = normalChromosomeRegionDepth / normalChromosomeRegionNum;
		System.err.println("normalChromosomeRegionNum:"+normalChromosomeRegionNum+"\trecalChrDepth:"+recalChrDepth);
		for(String chrName :  chrRd.keySet()) {
			System.err.println(chrName + "\t" + chrRd.get(chrName)/recalChrDepth);
		}
		
		Sex chrx = null;
		Sex chry = null;
		for(String chrName : sexChrsDepth.keySet()) {
			double formatDepth;
			
			if(isChr(chrName, "chrx")) {
				formatDepth = sexChrsDepth.get(chrName) / recalChrDepth;
				System.err.println("chrx formatDepth:" + formatDepth);
				if(!Double.isNaN(formatDepth)) {
					if(formatDepth > 0.75)
						chrx = Sex.F;
					else
						chrx = Sex.M;
				}
			}
			if(isChr(chrName, "chry")) {
				formatDepth = sexChrsDepth.get(chrName) / recalChrDepth;
				System.err.println("chry formatDepth:" + formatDepth);
				if(!Double.isNaN(formatDepth)) {
					if(formatDepth < 0.25)
						chry = Sex.F;
					else
						chry = Sex.M;
				}
			}
		}
		Sex gender = null;
		if(chrx == null && chry == null) {
			BAMQCException e = new BAMQCException.GenderInformationException("no gender chromosome info");
			e.printStackTrace();
			return Sex.unKown;
		}
		if(chrx != null && chry != null && chrx != chry) {
			BAMQCException e = new BAMQCException.GenderInformationException("gender is incompatible from chrX and chrY");
			e.printStackTrace();
			return Sex.unKown;
		}
		if(chrx != null && chry == null)
			gender = chrx;
		if(chrx == null && chry != null)
			gender = chry;
		if(chrx != null && chry != null && chrx == chry)
			gender = chrx;
		
		System.err.println("gender:" + gender.toString());
		return gender;
	}
	
	private void lowess() {
		//LoessInterpolator loess = new LoessInterpolator();
		Lowess lowess = new Lowess(0.1, 5);
		ArrayList<Double> depth = new ArrayList<Double>();
		ArrayList<Double> GC = new ArrayList<Double>();
		Map<Double, Double> gcDepth = new HashMap<>();
		
		for(String chrName : chrIndex.keySet()) {
			if(!isSexChromosome(chrName)) {
				Integer[] se = chrIndex.get(chrName);
				for(int i = se[0]; i <= se[1]; i++) {
					gcDepth.put(lowessData.get(i).getGCrate(), lowessData.get(i).getDepth());
				}
			}
		}
		
		Map<Double, Double> sGcDepth = new TreeMap<>(gcDepth);
		double basicGC = 0.0;
		ArrayList<Double> tmpDepth = new ArrayList<Double>();
		for(Map.Entry<Double, Double> data : sGcDepth.entrySet()) {
			if(basicGC == 0.0)
				basicGC = data.getKey();
			if(data.getKey() > basicGC + span) {
				double meanDepth = calMidVale(tmpDepth);
				GC.add(basicGC + span/2);
				depth.add(meanDepth);
				//System.err.println("add GC:" + basicGC + "-" + (basicGC + span) + "\tdepth:" + meanDepth);
				tmpDepth.clear();
				while(data.getKey() > basicGC + span) {
					basicGC = basicGC + span;
				}
			}
			//System.err.println("GC:" + data.getKey() + "\tdepth:" + data.getValue());
			tmpDepth.add(data.getValue());	
		}
		
		List<Double> recalDepth = lowess.lowess(GC, depth, depth.size());
		recaltable = new ArrayList<Double>();
		for( int i = 0; i < recalDepth.size(); i++) {
			//System.err.println("GC:" + GC.get(i) + "\tdepth:" + depth.get(i) + "\trecalDepth:" + recalDepth.get(i) + "\twight:" + totalMeanDepth / recalDepth.get(i));
			if(recalDepth.get(i) <= 0)
				continue;
			recaltable.add(totalMeanDepth / recalDepth.get(i));
			GCFinal.add(GC.get(i));
			//System.err.println("GC:" + GC.get(i) + "\tdepth:" + depth.get(i) + "\trecalDepth:" + recalDepth.get(i) + "\twight:" + totalMeanDepth / recalDepth.get(i));
		}
	}
	
	public double getRecalDepth(double GCValue, double depth) {
		int mid = GCFinal.size() / 2;
		//System.err.println("input:" + GCValue + "\t" + depth);
		for(int start = 0, end = GCFinal.size(); start <= end; mid = (start + end) / 2) {
			double beforeDis = 1, afterDis = 1, midDis = 1;
			//System.err.println("start:" + start + "\tend:" + end + "\tmid:" + mid);
			if(mid - 1 >= 0) 
				beforeDis = Math.abs(GCValue - GCFinal.get(mid - 1));
			if(mid + 1 < GCFinal.size())
				afterDis = Math.abs(GCValue - GCFinal.get(mid + 1));
			midDis = Math.abs(GCValue - GCFinal.get(mid));
			//System.err.println("before:" + beforeDis + "\tmid:" + midDis + "\tafter:" + afterDis);
			
			if(midDis > beforeDis)
				end = mid - 1;
			if(midDis > afterDis)
				start = mid + 1;
			if(midDis <= beforeDis && midDis <= afterDis) {
				//System.err.println("output:" + GCFinal.get(mid) +"\t" + recaltable.get(mid) * depth);
				return recaltable.get(mid) * depth;
			}
			//System.err.println("start2:" + start + "\tend2:" + end);
		}
		
		return -1;
	}
	
	private boolean isSexChromosome(String chrName) {
		String formatChrName = formatChrName(chrName);
		
		if(formatChrName.equals("chrx") || formatChrName.equals("chry"))
			return true;
		return false;
	}
	
	private boolean isChr(String chrName, String targetChrName) {
		String formatChrName = formatChrName(chrName);
		
		if(formatChrName.equals(targetChrName))
			return true;
		return false;
	}
	
	private String formatChrName(String chrName) {
		String formatChrName = chrName.toLowerCase();
		if(!formatChrName.startsWith("chr"))
			formatChrName = "chr"+formatChrName;
		return formatChrName;
	}
	
	public class LowessData {
		private double depth;
		private double GCrate;
		
		public LowessData(double depth, double GCrate) {
			this.depth = depth;
			this.GCrate = GCrate;
		}

		public double getDepth() {
			return depth;
		}

		public double getGCrate() {
			return GCrate;
		}
	}
	
	private double calMidVale(List<Double> list) {
		if(list == null || list.size() == 0)
			return 0;
		double valSum = 0;
		for(Double val : list)
			valSum += val;
		return valSum / list.size();
	}
	
	public class GCDepthComparator implements Comparator<Map.Entry<Double, Double>> {

		@Override
		public int compare(Entry<Double, Double> o1, Entry<Double, Double> o2) {
			double val1 = o1.getKey();
			double val2 = o2.getValue();
			
			if(val1 > val2)
				return 1;
			else if(val1 < val2)
				return -1;
			else
				return 0;
		}
		
	}
	
	public class LowessDataComparator implements Comparator<LowessData> {
		@Override
		public int compare(LowessData o1, LowessData o2) {
			if(o1.getDepth() > o2.getDepth())
				return 1;
			else if(o1.getDepth() < o2.getDepth())
				return -1;
			else
				return 0;
		}
	}
	
	public class LowessDataMapComparator implements Comparator<Entry<Integer, LowessData>> {
		@Override
		public int compare(Entry<Integer, LowessData> o1, Entry<Integer, LowessData> o2) {
			double val1 = o1.getValue().getGCrate();
			double val2 = o2.getValue().getGCrate();
			if(val1 > val2)
				return 1;
			else if(val1 < val2)
				return -1;
			else
				return 0;
		}
		
	}
}
