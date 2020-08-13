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
package org.bgi.flexlab.gaea.data.mapreduce.input.bed;

import org.bgi.flexlab.gaea.data.structure.region.TargetRegion;
import org.bgi.flexlab.gaea.util.FileIterator;

import java.io.IOException;

public class RegionHdfsParser extends TargetRegion{
	
	public void parseBedFileFromHDFS(String bedFilePath, boolean isWithFlank) throws IOException {
		FileIterator it = new FileIterator(bedFilePath);
		while(it.hasNext()) {
			parseBedRegion(it, isWithFlank);
		}
		it.close();
		//System.err.println("region size:" + regionSize);
	}
	
	protected void parseBedRegion(FileIterator it, boolean isWithFlank) {
		String line=it.next().toString().trim();
		String[] splitArray = line.split("\\s+");
		if(line.equals("") || line == null) {
			return;
		}
		
		boolean skipLine = parseBedFileLine(splitArray, this);
		if(skipLine) {
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
}
