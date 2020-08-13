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
package org.bgi.flexlab.gaea.data.mapreduce.input.fastq;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;

import java.io.IOException;

public class FastqRecordReader extends RecordReader<Text, Text> {
	public final static String READ_NAME_TYPE = "read.name.type";

	private FastqBasicReader reader = null;
	private Text key = null;
	private Text value = null;

	public FastqRecordReader() throws IOException {
	}

	@Override
	public void close() throws IOException {
		if (reader != null)
			reader.close();
	}

	@Override
	public float getProgress() throws IOException {
		return reader.getProgress();
	}

	@Override
	public Text getCurrentKey() throws IOException, InterruptedException {
		return key;
	}

	@Override
	public Text getCurrentValue() throws IOException, InterruptedException {
		return value;
	}

	@Override
	public void initialize(InputSplit split, TaskAttemptContext context)
			throws IOException, InterruptedException {
		Configuration configuration = context.getConfiguration();
		int readNameType = configuration.getInt(READ_NAME_TYPE, 0);
		byte[] recordDelimiter = null;
		if (configuration.get("textinputformat.record.delimiter") != null){
			recordDelimiter = configuration.get(
					"textinputformat.record.delimiter").getBytes();
		}
		if (readNameType == 0) {// read id format : reads_XX/1
			reader = new FastqForwardSlashReader(configuration,
					(FileSplit) split, recordDelimiter);
		} else if (readNameType == 1) {// read id format : reads_xx: 1:N:XX
										// reads_xx: 2:N:XX
			reader = new FastqSapceReader(configuration, (FileSplit) split,
					recordDelimiter);
		} else if (readNameType == 2) {// read id format : reads_xx
			reader = new FastqSpecialReader(configuration, (FileSplit) split,
					recordDelimiter);
		}
	}

	@Override
	public boolean nextKeyValue() throws IOException, InterruptedException {
		if(key == null)
			key = new Text();
		if(value == null)
			value = new Text();
		return reader.next(key, value);
	}
}
