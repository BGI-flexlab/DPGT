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
package org.bgi.flexlab.gaea.tools.mapreduce.bamqualitycontrol;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.commons.lang.math.RandomUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.SamRecordDatum;
import org.bgi.flexlab.gaea.util.SamRecordUtils;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class BamQualityControlMapper extends Mapper<LongWritable, SamRecordWritable, Text, Text>{
	/**
	 * FileHeader
	 */
	private SAMFileHeader mFileHeader=null;
	
	private String sampleName = "";
	
	private int unmappedReadsNum = 0;
	
	private int randomkey = RandomUtils.nextInt();
	
	private Text outK;
	
	private Text outV;
	
	private Map<String, Integer> rg2Index = new HashMap<String, Integer>();
	
	@Override
	public void setup(Context context) throws IOException {
		Configuration conf = context.getConfiguration();
		//header
		mFileHeader = SamHdfsFileHeader.getHeader(conf);
		assignIndexToReadGroup(mFileHeader);
	}

	@Override
	public void map(LongWritable key, SamRecordWritable value,Context context) throws IOException, InterruptedException {
		SamRecordDatum datum = new SamRecordDatum();
		String rgID = SamRecordUtils.getReadGroup(value.get());
		GaeaSamRecord record = new GaeaSamRecord(mFileHeader, value.get());
		sampleName = mFileHeader.getReadGroup(rgID).getSample();
		if(datum.parseSam(record.getSAMString())) {
			long winNum = -1;
			winNum = datum.getPosition() / BamQualityControl.WINDOW_SIZE;
			formatKeyValue(datum, rgID, winNum, true);
			context.write(outK, outV);
			if (winNum != (datum.getEnd() / BamQualityControl.WINDOW_SIZE)) {
				winNum++; 
				formatKeyValue(datum, rgID, winNum, false);
				context.write(outK, outV);
			}
		} else {
			if(unmappedReadsNum > 10000) {
				randomkey = RandomUtils.nextInt();
				unmappedReadsNum = 0;
			}
			context.write(new Text(formatKey("-1:-1", randomkey)), new Text("1"));
			unmappedReadsNum++;
		}
	}
	
	private void formatKeyValue(SamRecordDatum datum, String rgID, long winNum, boolean sameWin) {
		outK = new Text(formatKey(datum.getChrName(), winNum));
		outV = new Text(formatValue(datum, rgID, sameWin));
	}

	private String formatValue(SamRecordDatum datum, String rgID, boolean sameWin) {
		// TODO Auto-generated method stub
		StringBuilder vBuilder = new StringBuilder();
		vBuilder.append(datum.getFlag());
		vBuilder.append("\t");
		vBuilder.append(datum.getReadsSequence());
		vBuilder.append("\t");
		vBuilder.append(datum.getInsertSize());
		vBuilder.append("\t");
		vBuilder.append(datum.getPosition());
		vBuilder.append("\t");
		vBuilder.append(datum.getCigarString());
		vBuilder.append("\t");
		vBuilder.append(datum.getBestHitCount());
		vBuilder.append("\t");
		if(sameWin)
			vBuilder.append("0");
		else
			vBuilder.append("1");
		vBuilder.append("\t");
		vBuilder.append(datum.getMappingQual());
		vBuilder.append("\t");
		vBuilder.append(rg2Index.get(rgID));
		vBuilder.append("\t");
		vBuilder.append(datum.getQualityString());
		return vBuilder.toString();
	}
	
	private String formatKey(String chr, long winNum) {
		// TODO Auto-generated method stub
		StringBuilder kBuilder = new StringBuilder();
		kBuilder.append(sampleName);
		kBuilder.append(":");
		kBuilder.append(chr);
		kBuilder.append(":");
		kBuilder.append(winNum);
		return kBuilder.toString();
	}

	private void assignIndexToReadGroup(SAMFileHeader mFileHeader2) {
		// TODO Auto-generated method stub
		int rgIndex = 0;
		for(SAMReadGroupRecord rg : mFileHeader.getReadGroups()) {
			rg2Index.put(rg.getId(), rgIndex);
			rgIndex++;
		}
	}
}
