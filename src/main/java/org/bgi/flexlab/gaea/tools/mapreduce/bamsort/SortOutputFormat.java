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

import java.io.IOException;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.output.bam.GaeaBamOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.cram.GaeaCramOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;

public class SortOutputFormat extends
		FileOutputFormat<NullWritable, SamRecordWritable> {
	public static final String OUTPUT_NAME_PROP = "hadoopbam.sort.output.name",
			WRITE_HEADER_PROP = "hadoopbam.sort.output.write-header",
			OUTPUT_SAM_FORMAT_PROPERTY = "hadoopbam.anysam.output-format";

	private GaeaCramOutputFormat<NullWritable> cramOF = null;
	private GaeaBamOutputFormat<NullWritable> baseOF = null;

	private void initBaseOF(Configuration conf) {
		if (baseOF != null || cramOF != null)
			return;
		String format = conf.get(OUTPUT_SAM_FORMAT_PROPERTY);
		boolean multiSample = conf.getBoolean("multi.samples", true);
		if (format.toUpperCase().equals("CRAM")) {
			cramOF = new GaeaCramOutputFormat<>();
			cramOF.setHeader(multiSample);
		} else {
			baseOF = new GaeaBamOutputFormat<>();
			baseOF.setHeader(false);
//			baseOF.setHeader(conf.getBoolean(WRITE_HEADER_PROP,true));
		}
	}

	@Override
	public RecordWriter<NullWritable, SamRecordWritable> getRecordWriter(
			TaskAttemptContext context) throws IOException {
		initBaseOF(context.getConfiguration());

//		if (cramOF != null)
//			return cramOF.getRecordWriter(context,
//					getDefaultWorkFile(context, ""));

//		if (baseOF.getSAMHeader() == null) {
//			baseOF.setSAMHeader(SamHdfsFileHeader.getHeader(context.getConfiguration()));
//		}

		return baseOF.getRecordWriter(context, getDefaultWorkFile(context, ""));
	}

	@Override
	public void checkOutputSpecs(JobContext job) {
	}
}
