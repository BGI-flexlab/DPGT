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
package org.bgi.flexlab.gaea.data.structure.positioninformation.depth;

import org.bgi.flexlab.gaea.data.structure.positioninformation.CalculatePositionInforamtionInterface;
import org.bgi.flexlab.gaea.data.structure.positioninformation.CompoundInformation;
import org.bgi.flexlab.gaea.util.SamRecordDatum;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class PositionDepthCNV implements CalculatePositionInforamtionInterface<SamRecordDatum>{
	private ArrayList<PositionDepthCNVLane> cnvUsedDepth = new ArrayList<PositionDepthCNVLane>();
	private Map<Integer, Integer> rgIndex2depthIndex = new HashMap<Integer, Integer>();
	private int cnvDepthIndex = 0;
	private int laneSize;
	
	public PositionDepthCNV(int laneSize, int windowSize) {
		for(int i = 0; i < laneSize; i++) {
			PositionDepthCNVLane laneDepth = new PositionDepthCNVLane(windowSize);
			cnvUsedDepth.add(laneDepth);
		}
		this.laneSize = laneSize;
	}

	@Override
	public void add(CompoundInformation<SamRecordDatum> posInfo) {
		if(posInfo.eligiblePos() ) {
				//readInfo.getMappingQual() >= 10 && readInfo.getQualValue(readPosition) >= 15) {
			int depthIndex = 0;
			if(rgIndex2depthIndex.containsKey(posInfo.getRgIndex())) {
				depthIndex = rgIndex2depthIndex.get(posInfo.getRgIndex());
			} else {
				depthIndex = cnvDepthIndex;
				if(cnvDepthIndex >= laneSize){
					throw new RuntimeException("input data has more lane than BAM header!");
				}
				rgIndex2depthIndex.put(posInfo.getRgIndex(), cnvDepthIndex++);
			}
			PositionDepthCNVLane cnvLaneDepth = cnvUsedDepth.get(depthIndex);
			cnvLaneDepth.add(posInfo.distBetweenRefPosAndWinStart());
		}
	}
	
	public int getLaneSize() {
		return laneSize;
	}

	public PositionDepthCNVLane getLaneDepth(int index) {
		return cnvUsedDepth.get(index);
	}
}
