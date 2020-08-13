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
package org.bgi.flexlab.gaea.data.mapreduce.input.adaptor;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.CompressionCodecFactory;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.apache.hadoop.util.LineReader;

import java.io.IOException;

public class AdaptorRecordReader extends RecordReader<Text, Text>{
	protected static final Log LOG = LogFactory.getLog(AdaptorRecordReader.class.getName());
	protected CompressionCodecFactory compressionCodecs = null;
	private long start;
	protected long pos;
	protected long end;
	protected int maxLineLength;
	protected LineReader in;
	protected Text key = new Text();
	protected Text value = new Text();

	@Override
	public void initialize(InputSplit genericSplit,TaskAttemptContext context) throws IOException {
		FileSplit split = (FileSplit) genericSplit;
		System.out.println(split.toString());
		Configuration job = context.getConfiguration();
		System.err.println(split.getPath().toString());
		this.maxLineLength = job.getInt("mapred.linerecordreader.maxlength", Integer.MAX_VALUE);
		start = split.getStart();
	    end = start + split.getLength();
	    final Path file = split.getPath();
	    compressionCodecs = new CompressionCodecFactory(job);
	    final CompressionCodec codec = compressionCodecs.getCodec(file);

	    // open the file and seek to the start of the split
	    FileSystem fs = file.getFileSystem(job);
	    FSDataInputStream fileIn = fs.open(split.getPath());
	    boolean skipFirstLine = false;
	    if (codec != null) {
	      in = new LineReader(codec.createInputStream(fileIn), job);
	      end = Long.MAX_VALUE;
	    } else {
	      if (start != 0) {
	        skipFirstLine = true;
	        --start;
	        fileIn.seek(start);
	      }
	      in = new LineReader(fileIn, job);
	    }
	    if (skipFirstLine) {  // skip first line and re-establish "start".
	      start += in.readLine(new Text(), 0, (int)Math.min((long)Integer.MAX_VALUE, end - start));
	    }
	    this.pos = start;
	}
	
	@Override
	public boolean nextKeyValue() throws IOException {
		Text temp = new Text();
		while (pos < end) {
			int newSize = in.readLine(temp, maxLineLength, Math.max((int)Math.min(Integer.MAX_VALUE, end-pos), maxLineLength));
		    if (newSize == 0) {
		        return false;
		    }
		    pos += newSize;
		    if (newSize < maxLineLength) {
		    	String line = temp.toString();
		    	String[] splitLines = line.split("\t");
		    	if(splitLines[0].startsWith("#")) {
		    		continue;
		    	}
		    	String tempkey,tempvalue;
				int index = splitLines[0].lastIndexOf("/");
				if(index == -1){
					String[] str = splitLines[0].split(" ");
	        		tempvalue = str[1].substring(0, 1);
	        		tempkey = str[0]+"_1"+str[1].substring(1);
				}else{
					tempkey = splitLines[0].substring(0, index).trim();
	        		tempvalue = splitLines[0].substring(index + 1).trim();
				}
	        	
		    	key.set(tempkey);
				value.set(tempvalue);
		        return true;
		    }

		    // line too long. try again
		    LOG.info("Skipped line of size " + newSize + " at pos " + (pos - newSize));
		}
		    return false;
	}
	
	@Override
	public Text getCurrentKey() 
	{
		return key;
	}

	@Override
	public Text getCurrentValue()
	{
		return value;
	}

	@Override
	public void close() throws IOException {
		if (in != null) {
		   in.close(); 
		}
	}

	@Override
	public float getProgress() throws IOException, InterruptedException {
		if (start == end) {
			return 0.0f;
		} else {
			return Math.min(1.0f, (pos - start) / (float)(end - start));
		}
	}
}
