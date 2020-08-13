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

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.BlockLocation;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.CombineFileSplit;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;

import java.io.IOException;

public class GaeaCombineCramFileRecordReader extends
		RecordReader<LongWritable, SamRecordWritable> {
	protected GaeaCramRecordReader currentReader = null;
	protected CombineFileSplit split;
	protected int fileIndex;
	protected TaskAttemptContext context;
	
	public GaeaCombineCramFileRecordReader(InputSplit split, TaskAttemptContext context){
		this.split = (CombineFileSplit)split;
		this.context = context;
		this.fileIndex = 0;
		
		System.err.println("total number of file is "+this.split.getNumPaths());
		try {
			initializeNextRecordReader();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	@Override
	public void close() throws IOException {
		if (currentReader != null) {
			currentReader.close();
			currentReader = null;
		}
	}

	@Override
	public LongWritable getCurrentKey() throws IOException,
			InterruptedException {
		return currentReader.getCurrentKey();
	}

	@Override
	public SamRecordWritable getCurrentValue() throws IOException,
			InterruptedException {
		return currentReader.getCurrentValue();
	}

	@Override
	public float getProgress() throws IOException, InterruptedException {
		int subprogress = 0;
		if(fileIndex != 0)
			subprogress = fileIndex - 1;
		return (float) (subprogress) / split.getNumPaths();
	}

	@Override
	public void initialize(InputSplit split, TaskAttemptContext context)
			throws IOException, InterruptedException {
		if(null == currentReader){
			initializeNextRecordReader();
		}
	}

	@Override
	public boolean nextKeyValue() throws IOException, InterruptedException {
		while ((currentReader == null) || !currentReader.nextKeyValue()) {
			if (!initializeNextRecordReader()) {
				return false;
			}
		}
		return true;
	}

	protected boolean initializeNextRecordReader() throws IOException {
		if (currentReader != null) {
			currentReader.close();
			currentReader = null;
		}

		// if all chunks have been processed, nothing more to do.
		if (fileIndex == split.getNumPaths()) {
			return false;
		}

		// get a record reader for the fileIndex-th chunk
		try {
			Configuration conf = context.getConfiguration();
		    
			currentReader = new GaeaCramRecordReader();
			
			Path path = split.getPath(fileIndex);
			long length = split.getLength(fileIndex);
			FileSystem fs = path.getFileSystem(conf);
			FileStatus status = fs.getFileStatus(path);
			BlockLocation[] blkLocations = fs.getFileBlockLocations(status, 0, length);
			
			currentReader.initialize(new FileSplit(path, 0, length, blkLocations[0].getHosts()), context);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		fileIndex++;
		return true;
	}
}
