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
package org.bgi.flexlab.gaea.tools.mapreduce.bamsort;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;

public final class SortReducer
		extends
		Reducer<LongWritable, SamRecordWritable, NullWritable, SamRecordWritable> {
	private SortMultiOutputs<NullWritable, SamRecordWritable> mos;
	private SAMFileHeader header;
	private Map<String, String> formatSampleName = new HashMap<>();
	private String firstSample = null;
	private BamSortOptions options;

	@Override
	protected void setup(Context context) throws IOException {
		Configuration conf = context.getConfiguration();
		SAMFileHeader _header = SamHdfsFileHeader.getHeader(conf);;
		options = new BamSortOptions();
		options.getOptionsFromHadoopConf(conf);
		if(options.getRenames() != null)
			header = BamSortUtils.replaceSampleName(_header.clone(), options.getRenames());
		else
			header = _header;
		header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		for (SAMReadGroupRecord rg : header.getReadGroups()) {
			if (!formatSampleName.containsKey(rg.getSample()))
				formatSampleName.put(rg.getSample(),
						BamSortUtils.formatSampleName(rg.getSample()));
			if(firstSample == null)
				firstSample = BamSortUtils.formatSampleName(rg.getSample());
		}
		mos = new SortMultiOutputs<>(context);
	}

	public void replaceSampleNames(SAMFileHeader _header,String list) throws IOException{
		FileReader fr = new FileReader(list);
		BufferedReader br = new BufferedReader(fr);

		String line;
		HashMap<String,String> replaceList = new HashMap<>();
		while ((line = br.readLine()) != null) {
			String[] sampleNames = line.split("\t");
			replaceList.put(sampleNames[0], sampleNames[1]);
		}

		br.close();
		fr.close();
		header = BamSortUtils.replaceSampleName(_header.clone(), replaceList);
	}

	@Override
	protected void reduce(
			LongWritable ignored,
			Iterable<SamRecordWritable> records,
			Context ctx)
			throws IOException, InterruptedException {

		for (SamRecordWritable rec : records) {
			SAMRecord sam = rec.get();
			sam.setHeader(header);
			String sampleName = sam.getReadGroup().getSample();
			String fileName = formatSampleName.get(sampleName);
			if(!options.isMultiSample())
				fileName = firstSample;
			mos.write(fileName, NullWritable.get(), rec);

		}
	}

	@Override
	protected void cleanup(Context context) throws IOException,
			InterruptedException {
		mos.close();
	}
}
