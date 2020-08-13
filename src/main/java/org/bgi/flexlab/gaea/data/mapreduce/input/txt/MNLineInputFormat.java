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
package org.bgi.flexlab.gaea.data.mapreduce.input.txt;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.apache.hadoop.mapreduce.lib.input.NLineInputFormat;
import org.apache.hadoop.util.LineReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * MNLineInputFormat 继承 NLineInputFormat， 
 * 将输入数据按行平分到固定数目的mapper中
 */
public class MNLineInputFormat extends NLineInputFormat {
	public static final String MAPPER_NUM = "mapreduce.job.maps";
	
	/** 
	 * Logically splits the set of input files for the job, splits N lines
	 * of the input as one split.
	 * 
	 * @see FileInputFormat#getSplits(JobContext)
	 */
	@Override
	public List<InputSplit> getSplits(JobContext job) throws IOException {
		long start = System.currentTimeMillis();
		
		List<InputSplit> splits = new ArrayList<InputSplit>();
		int mapperNum = getMapperNum(job);
		long minNumLinesPerSplit = getMinNumLinesToSplit(job);
		long length = 0;
		List<FileStatus> fileList = listStatus(job);
		for (FileStatus file: fileList) {
		      length += file.getLen();
		}
		long meanLengthPerSplit = length/mapperNum;
		for (FileStatus status : fileList) {
			splits.addAll(getSplitsForFile(status,
					job.getConfiguration(), minNumLinesPerSplit, meanLengthPerSplit));
		}
		System.err.println("split耗时：" + (System.currentTimeMillis()-start)+"毫秒");		
		return splits;
	}
	
	public static List<FileSplit> getSplitsForFile(FileStatus status,
		      Configuration conf, long minNumLinesPerSplit, long meanLengthPerSplit) throws IOException {
		    List<FileSplit> splits = new ArrayList<FileSplit> ();
		    Path fileName = status.getPath();
		    if (status.isDirectory()) {
		      throw new IOException("Not a file: " + fileName);
		    }
		    FileSystem  fs = fileName.getFileSystem(conf);
		    LineReader lr = null;
		    try {
		      FSDataInputStream in  = fs.open(fileName);
		      lr = new LineReader(in, conf);
		      Text line = new Text();
		      int numLines = 0;
		      long begin = 0;
		      long length = 0;
		      int num = -1;
		      while ((num = lr.readLine(line)) > 0) {
		        numLines++;
		        length += num;
		        if (numLines >= minNumLinesPerSplit && length >= meanLengthPerSplit) {
		          splits.add(createFileSplit(fileName, begin, length));
		          begin += length;
		          length = 0;
		          numLines = 0;
		        }
		      }
		      if (numLines != 0) {
		        splits.add(createFileSplit(fileName, begin, length));
		      }
		    } finally {
		      if (lr != null) {
		        lr.close();
		      }
		    }
		    return splits;
		  }
	
	  public static List<FileSplit> getSplitsForFile(FileStatus status,
		      Configuration conf, long numLinesPerSplit) throws IOException {
		    List<FileSplit> splits = new ArrayList<FileSplit> ();
		    Path fileName = status.getPath();
		    if (status.isDirectory()) {
		      throw new IOException("Not a file: " + fileName);
		    }
		    FileSystem  fs = fileName.getFileSystem(conf);
		    LineReader lr = null;
		    try {
		      FSDataInputStream in  = fs.open(fileName);
		      lr = new LineReader(in, conf);
		      Text line = new Text();
		      int numLines = 0;
		      long begin = 0;
		      long length = 0;
		      int num = -1;
		      while ((num = lr.readLine(line)) > 0) {
		        numLines++;
		        length += num;
		        if (numLines == numLinesPerSplit) {
		          splits.add(createFileSplit(fileName, begin, length));
		          begin += length;
		          length = 0;
		          numLines = 0;
		        }
		      }
		      if (numLines != 0) {
		        splits.add(createFileSplit(fileName, begin, length));
		      }
		    } finally {
		      if (lr != null) {
		        lr.close();
		      }
		    }
		    return splits;
		  }
	
	 /**
	   * Set the number of mapper
	   * @param job the job to modify
	   * @param mapperNum the number of mapper
	   */
	  public static void setMapperNum(Job job, int mapperNum) {
	    job.getConfiguration().setInt(MAPPER_NUM, mapperNum);
	  }

	  /**
	   * Get the number of mapper
	   * @param job the job
	   * @return the number of mapper
	   */
	  public static int getMapperNum(JobContext job) {
	    return job.getConfiguration().getInt(MAPPER_NUM, 1);
	  }
	  
	  @Deprecated
	  public static void setNumLinesPerSplit(Job job, int numLines) {
		    job.getConfiguration().setInt(LINES_PER_MAP, numLines);
	  }
	  
	  @Deprecated
	  public static int getNumLinesPerSplit(JobContext job) {
		    return job.getConfiguration().getInt(LINES_PER_MAP, 1);
	  }
	  
	  /**
	   * Set the min number of lines per split
	   * @param job the job to modify
	   * @param numLines the number of lines per split
	   */
	  public static void setMinNumLinesToSplit(Job job, int numLines) {
	    job.getConfiguration().setInt(LINES_PER_MAP, numLines);
	  }
	  
	  /**
	   * Get the min number of lines per split
	   * @param job the job
	   * @return the number of lines per split
	   */
	  public static int getMinNumLinesToSplit(JobContext job) {
	    return job.getConfiguration().getInt(LINES_PER_MAP, 500);
	  }
}
