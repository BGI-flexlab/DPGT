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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2010 Aalto University 
 *
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.mapreduce.input.cram;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.Log.LogLevel;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;

import java.io.IOException;

public class GaeaCramInputFormat extends
		FileInputFormat<LongWritable, SamRecordWritable> {

	@Override
	public RecordReader<LongWritable, SamRecordWritable> createRecordReader(
			InputSplit inputSplit, TaskAttemptContext context)
			throws IOException, InterruptedException {
		Log.setGlobalLogLevel(LogLevel.ERROR);
		Configuration conf = context.getConfiguration();
		int type = conf.getInt("cram.inputformat.type", 0);
		CramInputFormatType inputType = CramInputFormatType.valueOf(type);
		
		RecordReader<LongWritable, SamRecordWritable> rr = null;
		if(inputType == CramInputFormatType.ALL)
			rr = new GaeaCramRecordReader();
		else if(inputType == CramInputFormatType.CHROMOSOME)
			rr = new GaeaCramChromosomeRecordReader();
		return rr;
	}
	
	public boolean isSplitable(JobContext job, Path path) {
		job.getConfiguration().setBoolean(GaeaCramRecordReader.CRAM_FILE_SPLITABLE,false);
		return false;
	}
}
