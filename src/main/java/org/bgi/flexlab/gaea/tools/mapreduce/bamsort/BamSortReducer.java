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
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.output.MultipleOutputs;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public final class BamSortReducer
		extends
		Reducer<LongWritable, SamRecordWritable, NullWritable, SamRecordWritable> {

	private SAMFileHeader samHeader;
	private MultipleOutputs<NullWritable,SamRecordWritable> mos;
//	private Map<Integer, String> sampleNames;

	@Override
	public void setup(Context context){
		Configuration conf = context.getConfiguration();
		samHeader = SamHdfsFileHeader.getHeader(conf);
		mos = new MultipleOutputs<>(context);

//		sampleNames = new HashMap<>();
//		List<SAMReadGroupRecord> list = samHeader.getReadGroups();
//		for(int i=0;i<list.size();i++)
//			sampleNames.put(i, list.get(i).getSample());
	}

	@Override
	protected void reduce(
			LongWritable key,
			Iterable<SamRecordWritable> records,
			Context ctx)
			throws IOException, InterruptedException {

		for (SamRecordWritable rec : records) {
			GaeaSamRecord sam = new GaeaSamRecord(samHeader, rec.get());
			SamRecordWritable w = new SamRecordWritable();
			w.set(sam);
			mos.write(NullWritable.get(), w, sam.getReadGroup().getSample());
		}
	}

	@Override
	protected void cleanup(Context context)
			throws IOException, InterruptedException {
		mos.close();
	}

}
