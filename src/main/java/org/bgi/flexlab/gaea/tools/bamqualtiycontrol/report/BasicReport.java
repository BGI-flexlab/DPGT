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

import org.apache.hadoop.mapreduce.Reducer.Context;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.BaseCounter;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.BaseType;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.ReadType;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.ReadsCounter;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.Tracker;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.Tracker.BaseTracker;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.Tracker.ReadsTracker;
import org.bgi.flexlab.gaea.util.SamRecordDatum;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;

public class BasicReport{
	
	private BaseTracker bTracker;
	
	private ReadsTracker rTracker;
	
	public BasicReport() {
		bTracker = new Tracker.BaseTracker();
		rTracker = new Tracker.ReadsTracker();
		register();
	}
	
	public boolean constructMapReport(SamRecordDatum datum, ReferenceShare genome, String chrName, Context context) {
		ChromosomeInformationShare chrInfo = genome.getChromosomeInfo(chrName);
		rTracker.setTrackerAttribute(ReadType.TOTALREADS);
		// 当位点坐标值+read长度大于染色体的长度时，则不处理该read，进入下一次循环
		if(datum.getEnd() >= chrInfo.getLength()) {
			if(context != null)
				context.getCounter("Exception", "read end pos more than chr end pos").increment(1);
			return false;
		}
		
		rTracker.setTrackerAttribute(ReadType.MAPPED);
		bTracker.setTrackerAttribute(BaseType.TOTALBASE.setCount(datum.getBaseCount()));
		
		if(datum.isUniqueAlignment()) {
			rTracker.setTrackerAttribute(ReadType.UNIQUE);
		}
			
		if ((datum.getFlag() & 0x400) != 0) {
			rTracker.setTrackerAttribute(ReadType.DUP);
		}
		
		if ((datum.getFlag() & 0x40) != 0 && (datum.getFlag() & 0x8) == 0) {
			rTracker.setTrackerAttribute(ReadType.PE);
		}
			
		String cigar = datum.getCigarString();
		if (cigar.contains("S") || cigar.contains("H")) {
			rTracker.setTrackerAttribute(ReadType.CLIPPED);
		}
		
		if (cigar.contains("D") || cigar.contains("I")) {
			rTracker.setTrackerAttribute(ReadType.INDEL);
		}
			
		if (isMismatch(datum, chrInfo)) {
			rTracker.setTrackerAttribute(ReadType.MISMATCHREADS);
		}
		return true;
	}
	
	@Override
	public String toString() {
		DecimalFormat df = new DecimalFormat("0.000");
		df.setRoundingMode(RoundingMode.HALF_UP);
		
		StringBuffer basicString = new StringBuffer();
		basicString.append("Basic Information:\n");
		basicString.append("Total Reads Number:\t");
		basicString.append(rTracker.getReadsCount(ReadType.TOTALREADS));
		basicString.append("\nTotal mapped bases:\t");
		basicString.append(bTracker.getProperty(BaseType.TOTALBASE));
		basicString.append("\nMapping Rate:\t");
		basicString.append(df.format(getRateOf(ReadType.MAPPED)));
		basicString.append("%\nUniq Mapping Rate:\t");
		if(getRateOf(ReadType.MAPPED) == getRateOf(ReadType.UNIQUE))
			basicString.append("NA");
		else
			basicString.append(df.format(getRateOf(ReadType.UNIQUE)));
		basicString.append("%\nPE Mapping Rate:\t");
		basicString.append(df.format(getRateOf(ReadType.PE)));
		basicString.append("%\nDuplication Rate:\t");
		basicString.append(df.format(getRateOf(ReadType.DUP)));
		basicString.append("%\nClipped Reads Rate:\t");
		basicString.append(df.format(getRateOf(ReadType.CLIPPED)));
		basicString.append("%\nMismatch Reads Rate:\t");
		basicString.append(df.format(getRateOf(ReadType.MISMATCHREADS)));
		basicString.append("%\nIndel Reads Rate:\t");
		basicString.append(df.format(getRateOf(ReadType.INDEL)));
		/*basicString.append("\nMismatch Rate in all mapped bases:\t");
		basicString.append(getMimatchRate());
		basicString.append("\nInsertion Rate in all mapped bases:\t");
		basicString.append(getInsertionRate());
		basicString.append("\nDeletion Rate in all mapped bases:\t");
		basicString.append(getDeletionRate());*/
		basicString.append("%\n");
		
		return basicString.toString();
	}
	
	public String toReducerString() {
		StringBuffer basicString = new StringBuffer();
		basicString.append("Basic Information:\n");
		for(Entry<String, ReadsCounter> counter : rTracker.getCounterMap().entrySet()) {
			basicString.append(counter.getKey());
			basicString.append(" ");
			basicString.append(counter.getValue().getReadsCount());
			basicString.append("\t");
		}
		for(Entry<String, BaseCounter> counter : bTracker.getCounterMap().entrySet()) {
			basicString.append(counter.getKey());
			basicString.append(" ");
			basicString.append(counter.getValue().getProperty());
			basicString.append("\t");
		}
		basicString.append("\n");
		
		return basicString.toString();
	}
	
	private boolean isMismatch(SamRecordDatum datum, ChromosomeInformationShare chrInfo) {
		boolean isMismatch = false;
		for (int i = (int) datum.getPosition(); i <= datum.getEnd(); i++) {
			int qpos = datum.getCigarState().resolveCigar(i, datum.getPosition());
			
			if (qpos < 0) {
				continue;
			}

			byte rbase = chrInfo.getBinaryBase(i);
			byte qbase = datum.getBinaryBase(qpos);

			if (rbase != qbase) {
				bTracker.setTrackerAttribute(BaseType.MISMATCH);
				isMismatch = true;
			}
		}
		datum.getCigarState().reSetCigarstate();
		return isMismatch;
	}
	
	public void parse(String line) {
		String[] splitArray = line.split("\t");
		for(String keyValue : splitArray)
			parseKeyValue(keyValue);
	}
	
	
	private void parseKeyValue(String keyValue) {
		String key = keyValue.split(" ")[0];
		String value = keyValue.split(" ")[1];
		ReadsCounter rCounter = null;
		BaseCounter bCounter = null;
		if((rCounter = rTracker.getCounterMap().get(key)) != null)
			rCounter.setReadsCount(Long.parseLong(value));
		else if((bCounter = bTracker.getCounterMap().get(key)) != null)
			bCounter.setBaseCount(Long.parseLong(value));
		else {
			throw new RuntimeException("Can not idenity counter with name " + key);
		}
			
	}
	
	public void register() {
		rTracker.register(createReadsCounters());
		bTracker.register(createBaseCounters());
	}
	
	public BaseTracker getBaseTracker() {
		return bTracker;
	}
	
	public ReadsTracker getReadsTracker() {
		return rTracker;
	}
	
	public double getRateOf(ReadType type) {
		long totalReads = rTracker.getReadsCount(ReadType.TOTALREADS);
		long target = rTracker.getReadsCount(type);
		return (100 * (target/(double)totalReads));
	}
	
	public List<BaseCounter> createBaseCounters() {
		List<BaseCounter> counters = new ArrayList<>();
		counters.add(new BaseCounter(BaseType.TOTALBASE));
		counters.add(new BaseCounter(BaseType.MISMATCH));
		return counters;
	}
	
	public List<ReadsCounter> createReadsCounters() {
		List<ReadsCounter> counters = new ArrayList<>();
		for(ReadType type : ReadType.values())
			counters.add(new ReadsCounter(type));
		return counters;
	}

}
