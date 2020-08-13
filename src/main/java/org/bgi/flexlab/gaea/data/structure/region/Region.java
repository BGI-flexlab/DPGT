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

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class Region {
	
	protected final static String INDEX_SEPARATOR = ".";
	/**
	 * 数据结构，方便查找:winID->[start,end],[start,end]...
	 */
	protected Map<String, ArrayList<Integer[]>> index;
	/**
	 * 染色体标志，表示此条染色体已经全在区域内
	 */
	protected ArrayList<String> chrs;
	/**
	 * 区域大小
	 */
	protected int regionSize;
	
	protected int start;
	
	protected int end;
	
	protected String chrName;
	
	public Region() {
		chrs = new ArrayList<String>();
		index = new ConcurrentHashMap<String, ArrayList<Integer[]>>();
		regionSize = 0;
	}
	
	public boolean isPositionInRegion(String chrName, long position) {
		if(chrs.contains(chrName)) {
			return true;
		}
		
		String id = formartID(position, chrName);
		
		if(index.containsKey(id)) {
			//System.out.println(id.toString());
			int regionNum = index.get(id).size();
			for(int i = 0; i < regionNum; i++) {
				if(position >= index.get(id).get(i)[0] && 
						position <= index.get(id).get(i)[1]){
					return true;
				}
			}
		}
		return false;
	}

	public boolean isSamRecordInRegion(SAMRecord samrecord) {
		return isReadInRegion(samrecord.getReferenceName(), samrecord.getStart() - 1, samrecord.getEnd() - 1);
	}
	
	public boolean isReadInRegion(String chrName, long start, long end) {
		if(isPositionInRegion(chrName, start) || isPositionInRegion(chrName, end))
			return true;
		else {
			for(long i = start + 1; i < end; i++) {
				if(isPositionInRegion(chrName, i))
					return true;
			}
		}
		return false;
	}
	
	public String formartID(long position, String chrName) {
		StringBuffer id = new StringBuffer();
		int winid = (int) position / 1000;
		id.append(chrName);
		id.append(INDEX_SEPARATOR);
		id.append(winid);
		return id.toString();
	}
	
	public Map<String, ArrayList<Integer[]>> getIndex(){
		return index;
	}
	
	public void updateRegionSize(int regionSize) {
		this.regionSize += regionSize;
	}
	
	public int getRegionSize() {
		return regionSize;
	}
	
	public void setChrName(String chrName) {
		this.chrName = chrName;
	}
	
	public String getChrName() {
		return chrName;
	}
	
	public void setStart(int start) {
		this.start = start;
	}
	
	public int getStart() {
		return start;
	}
	
	public void setEnd(int end) {
		this.end = end;
	}
	
	public int getEnd() {
		return end;
	}

}
