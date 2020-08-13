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

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.CompressionCodecFactory;
import org.apache.hadoop.io.compress.CompressionOutputStream;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegion;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegion.Regiondata;
import org.bgi.flexlab.gaea.tools.mapreduce.bamqualitycontrol.BamQualityControlOptions;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class CNVDepthReport {
	private static Map<String, Map<Integer, Integer>> bedIndex = null;
	private LaneDepth2[] depths;
	private int currentLaneID;
	private String currentChrName;
	
	public CNVDepthReport(int laneSize, SingleRegion region) {
		if(bedIndex == null) {
			bedIndex = new HashMap<String, Map<Integer,Integer>>();

			for(Regiondata data : region.getRegions()) {
				Map<Integer,Integer> bedRegionIndex;
				int index;
				int start = data.getStart();
				int end = data.getEnd();
				if(bedIndex.containsKey(data.getChrName())) {
					bedRegionIndex = bedIndex.get(data.getChrName());
					index = bedIndex.get(data.getChrName()).size() * 2;
					for(int pos = start; pos <= end; pos++) {
						if(!bedRegionIndex.containsKey(pos)) {
							bedRegionIndex.put(pos, index);
							index += 2;
						}
					}
				} else {
					bedRegionIndex = new HashMap<Integer, Integer>();
					index = 0;
					for(int pos = start; pos <= end; pos++) {
						if(!bedRegionIndex.containsKey(pos)) {
							bedRegionIndex.put(pos, index);
							index += 2;
						}
					}
					bedIndex.put(data.getChrName(), bedRegionIndex);
				}
			}
		}
		
		depths = new LaneDepth2[laneSize];
		for(int i = 0; i < laneSize; i++) {
			depths[i] = new LaneDepth2();
		}
	}
	
	/**
	 * result add
	 * @param line
	 */
	public void add(String line) {
		if(line.startsWith(">LaneID")) {
			currentLaneID = Integer.parseInt(line.split(":")[1]);
			return;
		}
		if(line.startsWith(">ChrName")) {
			currentChrName = line.split(":")[1];
			return;
		}
		String[] lineSplits = line.split("\t");
		int position = Integer.parseInt(lineSplits[0]);
		short depth = Short.parseShort(lineSplits[1]);
		depths[currentLaneID].add(currentChrName, position, depth);
	}
	
	public void toReport(BamQualityControlOptions options, FileSystem fs, Configuration conf, String sampleName) throws IOException {
		for(int i = 0; i < depths.length; i++) {
			Map<String, WrappedIntArray> sampleDepth = depths[i].laneDepth;
			for(String chrName : depths[i].laneDepth.keySet()) {
				StringBuffer cnvDepthFilePath = new StringBuffer();
				cnvDepthFilePath.append(options.getOutputPath());
				cnvDepthFilePath.append("/");
				cnvDepthFilePath.append("cnvDepth");
				cnvDepthFilePath.append("/");
				cnvDepthFilePath.append(sampleName);
				cnvDepthFilePath.append("-lane");
				cnvDepthFilePath.append(i);
				cnvDepthFilePath.append("/");
				cnvDepthFilePath.append(chrName);
				cnvDepthFilePath.append(".dep.gz");
				Path cnvDepthPath = new Path(cnvDepthFilePath.toString());
				FSDataOutputStream cnvDepthStream = fs.create(cnvDepthPath);
				CompressionCodecFactory codecFactory = new CompressionCodecFactory(conf);
		        CompressionCodec codec = codecFactory.getCodec(cnvDepthPath);
		        CompressionOutputStream compressedOutput = codec.createOutputStream(cnvDepthStream);
		        //ChrLaneDepth laneChrDepths = depths[i].laneDepth.get(chrName);
		        //Map<Integer, Integer> depthLanePos = laneChrDepths.depth;
		        int[] depth = sampleDepth.get(chrName).getArray();
		        StringBuilder sb = new StringBuilder();
		        for(int j = 0; j < depth.length; j += 2) {
		        	sb.append(chrName);
		        	sb.append("\t");
					sb.append(depth[j] + 1);
					sb.append("\t");
					sb.append(depth[j + 1]);
					sb.append("\n");
				}
		        compressedOutput.write(sb.toString().getBytes());
		        compressedOutput.close();
		        cnvDepthStream.close();
			}
		}
	}
	
	/**
	 * reducer add
	 * @param chrName
	 * @param position
	 * @param laneDepth
	 */
	public void add(String chrName, int position, int[] laneDepth) {
		if(laneDepth.length != depths.length)
			throw new RuntimeException("input lane size is wrong");
		for(int i = 0; i < laneDepth.length; i++) {
			depths[i].add(chrName, position, laneDepth[i]);
		}
	}
	
	/**
	 * to reducer output
	 * @return
	 */
	public String toReducerString() {
		StringBuilder sb = new StringBuilder();
		sb.append("CNV Depth:\n");
		for(int i = 0; i < depths.length; i++) {
			sb.append(">LaneID:");
			sb.append(i);
			sb.append("\n");
			sb.append(depths[i].toString());
		}
		sb.append("CNV Depth\n");
		return sb.toString();
	}
	
	public class LaneDepth2 {
		private Map<String, WrappedIntArray> laneDepth = new HashMap<String, WrappedIntArray>();
		
		public LaneDepth2() {
			for(String chrName : bedIndex.keySet()) {
				WrappedIntArray depth = new WrappedIntArray(bedIndex.get(chrName).size() * 2);
				
				laneDepth.put(chrName, depth);
			}
		}
		
		public void add(String chrName, int position, int depth) {
			int index = bedIndex.get(chrName).get(position);
			if(index < 0) {
				throw new RuntimeException("index < 0 when index position of CNV depth.");
			}
			if(index > laneDepth.get(chrName).size()) {
				System.err.println("chr:" + chrName + "\tposition:" + position + " OI");
			}
			laneDepth.get(chrName).set(index, position);
			laneDepth.get(chrName).set(index + 1, depth);
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder();
			for(String chrName : laneDepth.keySet()) {
				sb.append(">ChrName:");
				sb.append(chrName);
				sb.append("\n");
				//sb.append(laneDepth.get(chrName).toString());
				int[] depths = laneDepth.get(chrName).getArray();
				for(int index = 0; index < depths.length; index += 2) {
					if(depths[index + 1] != 0) {
						sb.append(depths[index]);
						sb.append("\t");
						sb.append(depths[index + 1]);
						sb.append("\n");
					}
				}
			}
			return sb.toString();
		}
	}
	
	public class WrappedIntArray {
		int[] array;
		
		public WrappedIntArray(int size) {
			array = new int[size];
			Arrays.fill(array, 0);
		}
		
		public int[] getArray() {
			return array;
		}
		
		public int size() {
			return array.length;
		}
		
		public void set(int index, int value) {
			array[index] = value;
		}
		
		public void add(int index, int value) {
			array[index] += value;
		}
	}
}
