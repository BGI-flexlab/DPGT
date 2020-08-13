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

import org.bgi.flexlab.gaea.util.ChromosomeUtils;
import org.bgi.flexlab.gaea.util.FileIterator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class TargetRegion extends Region {

	
	/**
	 * 每条染色体在区域中的大小
	 */
	private Map<String, Integer> chrsSize;

	protected FlankRegion flankRegion;
	
	public TargetRegion() {
		super();
		chrsSize = new ConcurrentHashMap<String, Integer>();
		flankRegion = new FlankRegion();
	}
	
	public void addRegionIndex(String chrName, int start, int end, boolean isAddFlank) {
		for(int i = (int) start/1000; i <= (int) end/1000; i ++) {
			int winStart, winEnd;
			
			winStart = initWinStart(start, i);
			winEnd = initWinEnd(end, i);
		
			StringBuffer id = new StringBuffer();
			id.append(chrName);
			id.append(INDEX_SEPARATOR);
			id.append(i);
			
			Integer[] winse = new Integer [2];
			winse[0] = winStart;
			winse[1] = winEnd;
			if(!isAddFlank) {
				addRegionIndex(this.index, id, winse);
			} else {
				addRegionIndex(flankRegion.getIndex(), id, winse);
			}
		}
	}
	
	private void addRegionIndex(Map<String, ArrayList<Integer[]>> index, StringBuffer id,
			Integer[] winse) {
		if(!index.containsKey(id.toString())) {
			ArrayList<Integer[]> wins = new ArrayList<Integer[]>();
			wins.add(winse);
			index.put(id.toString(), wins);
		}
		else {
			index.get(id.toString()).add(winse);
		}
	}
	
	private int initWinStart(int start, int i) {
		int winStart;
		if(i == start/1000) {
			winStart = start;
		} else {
			winStart = i * 1000;
		}
		return winStart;
	}
	
	private int initWinEnd(int end, int i) {
		int winEnd;
		if(i == end/1000) {
			winEnd = end;
		} else {
			winEnd = i * 1000 + 999;
		}
		return winEnd;
	}
	
	/**
	 * 解析命令行区域参数，chr:start-end
	 * @param region
	 */
	public void parseRegion(String region, boolean isWithFlank) {
		String[] splitArray = region.split(":");
		if(splitArray.length == 1) {
			System.out.println("region value is the whole chromosome.");
			chrs.add(splitArray[0]);
			return;
		}
		this.setChrName(splitArray[0]);

		String[] splitArray2 = splitArray[1].split("-");
		if(splitArray2.length == 1) {
			System.out.println("you must set an end point value.");
			return;
		}
		this.setStart(Integer.parseInt(splitArray2[0]) - 1);
		this.setEnd(Integer.parseInt(splitArray2[1]) - 1);
		flankRegion.setRegionSize(calRegionSize(start, end));
		
		addRegionIndex(chrName, start, end, false);
		if(isWithFlank) {
			addRegionIndex(chrName, start - flankRegion.getExtendSize() > 0 ?
					start - flankRegion.getExtendSize() : 0, 
					end + flankRegion.getExtendSize(), true);
		}
	}
	
	private void positionInRegion(String chrName,long position,ArrayList<Integer[]> list,HashSet<String> starts){
		String id = formartID(position, chrName);
		
		if(index.containsKey(id)) {
			//System.out.println(id.toString());
			int regionNum = index.get(id).size();
			int _start,_end;
			for(int i = 0; i < regionNum; i++) {
				_start = index.get(id).get(i)[0];
				_end = index.get(id).get(i)[1];
				if(_start > position)break;
				if(position >= _start && position <= _end){
					String k = _start+","+_end;
					if(starts.contains(k))continue;
					else{
						Integer[] tmp = new Integer[2];
						tmp[0] = _start;
						tmp[1] = _end;
						list.add(tmp);
						starts.add(k);
					}
				}
			}
		}
	}
	
	public ArrayList<Integer[]> readInRegion(String chrName,long start,long end){
		ArrayList<Integer[]> list = new ArrayList<Integer[]>();
		HashSet<String> starts = new HashSet<String>();
		
		Integer[] wins = new Integer[2];
		if(chrs.contains(chrName)) {
			wins[0] = wins[1] = -1;
			list.add(wins);
			return list;
		}
		
		positionInRegion(chrName,start,list,starts);
		positionInRegion(chrName,end,list,starts);

		starts.clear();
		return list;
	}
	
	public int getChrSize(String chrName) {
		//SSystem.out.println(chrName);
		chrName = ChromosomeUtils.formatChrName(chrName);
		//System.out.println(chrName);
		if(chrsSize.containsKey(chrName)) {
			return chrsSize.get(chrName);
		} else {
			return 0;
		}
	}
	
	protected void addChrSize(String chrName, int size) {
		chrName = ChromosomeUtils.formatChrName(chrName);
		if(chrsSize.containsKey(chrName)) {
			chrsSize.put(chrName, chrsSize.get(chrName) + size);
		} else {
			chrsSize.put(chrName, size);
		}
	}
	
	public void parseBedFileFromHDFS(String bedFilePath, boolean isWithFlank) throws IOException {
		FileIterator it = new FileIterator(bedFilePath);
		while(it.hasNext()) {
			parseBedRegion(it, isWithFlank);
		}

		it.close();
	}
	
	protected void parseBedRegion(FileIterator it, boolean isWithFlank) {
		String line=it.next().toString().trim();
		String[] splitArray = line.split("\\s+");
		if(line.equals("") || line == null) {
			return;
		}
		
		boolean skipLine = parseBedFileLine(splitArray, this);
		if(skipLine) {
			System.out.println("skip line");
			return;
		} else {
			addChrSize(chrName, calRegionSize(start, end));
		}
		
		if(isWithFlank) {
			while(it.hasNext()) {
				line=it.next().toString().trim();
				splitArray = line.split("\t");
				skipLine = parseBedFileLine(splitArray, flankRegion);
				if(skipLine) {
					continue;
				} else {
					addChrSize(flankRegion.getChrName(), calRegionSize(flankRegion.getStart(), flankRegion.getEnd()));
				}
				
				processWindow();
			}
			//处理最后一个窗口
			processWindow(start, end, chrName);
		}
	}
	
	@Deprecated
	public void parseBedFile(String bedFilePath, boolean isWithFlank) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(bedFilePath));
		String line = "";
		String[] splitArray;

		while((line = br.readLine()) != null) {
			splitArray = line.split("\t");
			boolean skipLine = parseBedFileLine(splitArray, this);
			if(skipLine) {
				continue;
			}
			if(isWithFlank) {
				while((line = br.readLine()) != null) {
					splitArray = line.split("\t");
					skipLine = parseBedFileLine(splitArray, flankRegion);
					if(skipLine) {
						continue;
					}
					processWindow();
				}
				//处理最后一个窗口
				processWindow(start, end, chrName);
			}
		}
		
		br.close();
	}
	
	protected boolean parseBedFileLine(String[] splitArray, Region region){
		if(splitArray.length == 1) {
			System.out.println("region value is the whole chromosome:" + splitArray[0]);
			chrs.add(splitArray[0]);
			return true;
		}
		if(splitArray.length == 2) {
			System.out.println("you must set an end point value.");
			return true;
		}
		region.setChrName(splitArray[0]);
		region.setStart(Integer.parseInt(splitArray[1]));
		region.setEnd(Integer.parseInt(splitArray[2]) - 1);
		this.updateRegionSize(calRegionSize(region.getStart(), region.getEnd()));
		addRegionIndex(region.getChrName(), region.getStart(), region.getEnd(), false);
		return false;
	}
	
	protected void processWindow() {
		if(notSorted()) {
			throw new RuntimeException("BED file is not sorted:" + flankRegion.getStart() + "<" + start);
		}

		if(combineWin()) {//窗口1和2合并
			this.setEnd(flankRegion.getEnd());
		} 
		
		if(startProcess()) {
			processWindow(start, end, chrName);
			this.setChrName(flankRegion.getChrName());
			this.setStart(flankRegion.getStart());
			this.setEnd(flankRegion.getEnd());
		}
	}
	
	protected void processWindow(int start, int end, String chrName){
		int flankStart = start - flankRegion.getExtendSize() > 0 ? start - flankRegion.getExtendSize() : 0;
		int flankEnd =  end + flankRegion.getExtendSize();
		flankRegion.updateRegionSize(calRegionSize(flankStart, flankEnd));
		addRegionIndex(chrName, flankStart, flankEnd, true);
	}

	
	protected int calRegionSize(int start, int end) {
		return (end - start + 1);
	}
	
	private boolean notSorted(){
		return this.chrName.equals(flankRegion.getChrName()) && flankRegion.getStart() < this.start;
	}
	
	private boolean combineWin() {
		return this.end + flankRegion.getExtendSize() >= flankRegion.getStart() - flankRegion.getExtendSize()
				&& this.chrName.equals(flankRegion.getChrName());
	}
	
	private boolean startProcess() {
		return this.end + flankRegion.getExtendSize() < flankRegion.getStart() - flankRegion.getExtendSize() 
				|| !this.chrName.equals(flankRegion.getChrName());
	}
	
	public  Map<String, ArrayList<Integer[]>> getFlankIndex() {
		return flankRegion.getIndex();
	}
	
	public FlankRegion getFlankRegion() {
		return flankRegion;
	}
	
	public int getRegionFlankSize() {
		return flankRegion.getRegionSize();
	}
	
	public boolean isReadInFlank(String chrName, long start, long end) {
		return flankRegion.isReadInRegion(chrName, start, end);
	}
	
	public boolean isPositionInFlank(String chrName, long position) {
		return flankRegion.isPositionInRegion(chrName, position);
	}
}
