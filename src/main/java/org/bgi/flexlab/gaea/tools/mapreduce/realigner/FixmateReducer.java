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
package org.bgi.flexlab.gaea.tools.mapreduce.realigner;

import htsjdk.samtools.SAMFileHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.GaeaSamPairUtil;

import java.io.IOException;

public class FixmateReducer extends
		Reducer<Text, SamRecordWritable, NullWritable, SamRecordWritable> {
	private SAMFileHeader header = null;
	private SamRecordWritable valueout = null;

	protected void setup(Context context) throws IOException,
			InterruptedException {
		Configuration conf = context.getConfiguration();
		header = SamHdfsFileHeader.getHeader(conf);
		valueout = new SamRecordWritable();
	}

	@Override
	public void reduce(Text key, Iterable<SamRecordWritable> values,
			Context context) throws IOException, InterruptedException {

		GaeaSamRecord[] reads = new GaeaSamRecord[2];
		reads[0] = null;
		reads[1] = null;

		for (SamRecordWritable read : values) {
			GaeaSamRecord sam = new GaeaSamRecord(header, read.get());
			int type = (((sam.getFlags() >> 6) & 1) == 1) ? 0 : 1;

			reads[type] = reads[type] == null ? sam : getBetterRead(
					reads[type], sam);
		}

		if (reads[0] == null || reads[1] == null) {
			for (int i = 0; i < 2; i++) {
				if (reads[i] != null)
					writeRead(reads[i], context);
			}
		} else {
			if (!pairedReads(reads[0], reads[1]))
				GaeaSamPairUtil.setMateInfo(reads[0], reads[1], header);

			writeRead(reads[0], context);
			writeRead(reads[1], context);
		}

		reads[0] = null;
		reads[1] = null;
	}

	private GaeaSamRecord getBetterRead(GaeaSamRecord last, GaeaSamRecord curr) {
		if (last.getMappingQuality() > curr.getMappingQuality())
			return last;
		else if (last.getMappingQuality() < curr.getMappingQuality())
			return curr;
		if (last.getFlags() < curr.getFlags())
			return curr;
		return last;
	}

	public boolean pairedReads(GaeaSamRecord rec1, GaeaSamRecord rec2) {
		if (rec1.getMateReferenceIndex() != rec2.getReferenceIndex())
			return false;
		if (rec1.getMateAlignmentStart() != rec2.getAlignmentStart())
			return false;
		if (rec2.getMateReferenceIndex() != rec1.getReferenceIndex())
			return false;
		if (rec2.getMateAlignmentStart() != rec1.getAlignmentStart())
			return false;
		return true;
	}

	private void writeRead(GaeaSamRecord read, Context context) {
		valueout.set(read);
		try {
			context.write(NullWritable.get(), valueout);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
}
