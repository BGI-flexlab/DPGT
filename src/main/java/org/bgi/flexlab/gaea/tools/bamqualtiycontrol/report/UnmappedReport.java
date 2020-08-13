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

import htsjdk.samtools.SAMRecordIterator;
import org.apache.hadoop.io.Text;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.ReadType;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class UnmappedReport {
	
	private Map<String, ArrayList<Long>> unmappedSites = new ConcurrentHashMap<String, ArrayList<Long>>();
	
	private long unmappedEnd = 0, unmappedStart = 0;
	
	private boolean firstUnmappedSite = true;
	
	public String toReducerString() {
		StringBuilder info = new StringBuilder();
		for(String key : unmappedSites.keySet()) {
			ArrayList<Long> unmaped = unmappedSites.get(key);
			for(int i = 0; i < unmaped.size(); i += 2) {
				info.append(unmaped.get(i));
				info.append("\t");
				info.append(unmaped.get(i + 1));
				info.append("\n");
			}
		}
		return info.toString();
	}
	
	public Map<String, ArrayList<Long>> getUnmappedSites() {
		return unmappedSites;
	}
	
	public ArrayList<Long> getUnmappedSites(String chrName) {
		if(chrName == null || chrName == "") {
			return null;
		}
		
		if(!unmappedSites.containsKey(chrName)) {
			ArrayList<Long> sites = new ArrayList<Long>();
			unmappedSites.put(chrName, sites);
		}
		
		return unmappedSites.get(chrName);
	}
	
	public boolean constructMapReport(long winNum, String chrName, Iterable<Text> values, BasicReport basicReport) {
		if(winNum < 0 || chrName.equals("-1")) {//unmapped
			Iterator<Text> vals = values.iterator();
			while (vals.hasNext()) {
				basicReport.getReadsTracker().setTrackerAttribute(ReadType.TOTALREADS);
				vals.next();
			}
			return true;
		}
		return false;
	}

	public boolean constructMapReport(long winNum, String chrName, SAMRecordIterator values, BasicReport basicReport) {
		if(winNum < 0 || chrName.equals("-1")) {//unmapped
			while (values.hasNext()) {
				basicReport.getReadsTracker().setTrackerAttribute(ReadType.TOTALREADS);
				values.next();
			}
			return true;
		}
		return false;
	}
	
	public void finalize(ArrayList<Long> unmappedSites) {
		if(unmappedEnd != 0 && unmappedStart != 0) {
			unmappedSites.add(unmappedStart);
			unmappedSites.add(unmappedEnd + 1);
		}
	}
	
	public void updateUnmappedSites(long pos, ArrayList<Long> unmappedSites) {
		if(firstUnmappedSite) {
			unmappedEnd = unmappedStart = pos;
			firstUnmappedSite = false;
		} else {
			if(pos == unmappedEnd + 1) {
				unmappedEnd++;
			} else {
				unmappedSites.add(unmappedStart);
				unmappedSites.add(unmappedEnd + 1);
				unmappedEnd = unmappedStart = pos;
			}
		}
	}

}
