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
package org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter;

import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.BaseType;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.Depth;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.DepthType;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.Interval;

public class BaseCounter {

	private long baseCount;
	
	private long baseWithoutPCRdupCount;
	
	private long totalDepth;
		
	private long totalDepthWithouPCR;
	
	private Interval region;
	
	private Depth depth;
	
	private DepthType dType;
	
	private BaseType bType;
	
	private CounterProperty[] properties;
	
	public BaseCounter(CounterProperty... counterProperties) {
//		if counter properties has 3 elements, this counter is working for the region report
		properties = counterProperties;
		if(properties.length == 3) {
			region = (Interval)properties[0];
			depth = (Depth)properties[1];
			dType = (DepthType)properties[2];
			bType = null;
		}
//		if counter properties has 1 element, this counter is working for the basic report
		if(properties.length == 1) {
			bType = (BaseType) properties[0];
		}
	}
	
	public void update(BaseType bType) {
		if(this.bType == bType) {
			this.baseCount += bType.getCount();
		}
	}
	
	public void update(Interval region, DepthType depth, DepthType noPCRdepth) {
		// TODO Auto-generated method stub
		if(this.region == region) {
			if(this.depth != Depth.TOTALDEPTH) {
				if(this.depth == Depth.ABOVE_ZREO && depth.getDepth() == 0) {
					baseCount++;
					baseWithoutPCRdupCount++;
				} else if(this.dType == DepthType.NORMAL && depth.getDepth() >= this.depth.getDepth()) {
					baseCount++;
				} else if(this.dType == DepthType.WITHOUT_PCR && noPCRdepth.getDepth() >= this.depth.getDepth()) {
					baseWithoutPCRdupCount++;
				}
			} else {
				switch (this.dType) {
				case NORMAL:
					totalDepth += depth.getDepth();
					break;
				case WITHOUT_PCR:
					totalDepthWithouPCR += noPCRdepth.getDepth();
				default:
					break;
				}
					
			 }
		}
	}

	public long getBaseCount() {
		return baseCount;
	}
	
	public void setBaseCount(long baseCount) {
		this.baseCount += baseCount;
	}
	
	public long getBaseWithoutPCRDupCount() {
		return baseWithoutPCRdupCount;
	}

	public void setBaseWithoutPCRDupCount(long baseWithoutPCRdupCount) {
		this.baseWithoutPCRdupCount += baseWithoutPCRdupCount;
	}
	
	public long getTotalDepth() {
		// TODO Auto-generated method stub
		return totalDepth;
	}
	
	public void setTotalDepth(long totalDepth) {
		this.totalDepth += totalDepth;
	}
	
	public void setTotalDepthWithoutPCRDup(long totalDepth) {
		this.totalDepthWithouPCR += totalDepth;
	}
	
	public long getProperty() {
		long result = 0;
		if(properties.length == 1) {
			result = baseCount;
			return result;
		} else {
			 switch (depth) {
			 case TOTALDEPTH:
				 switch (dType) {
				case NORMAL:
					result = totalDepth;
					break;
				case WITHOUT_PCR:
					result = totalDepthWithouPCR;
					break;
				}
				 break;
			 default:
				 switch (dType) {
				 case NORMAL:
				 	result = baseCount;
				 	break;
				 case WITHOUT_PCR:
					result = baseWithoutPCRdupCount;
				 	break;
				 }
			  }
			  return result;
		}
	}
	
	public String formatKey() {
		String key = "";
		for(CounterProperty property : properties)
			key += property.toString();
		return key;
	}

}
