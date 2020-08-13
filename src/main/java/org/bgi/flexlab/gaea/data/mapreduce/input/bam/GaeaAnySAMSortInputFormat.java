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
package org.bgi.flexlab.gaea.data.mapreduce.input.bam;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.bgi.flexlab.gaea.data.mapreduce.input.sam.GaeaSamInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.seqdoop.hadoop_bam.FileVirtualSplit;
import org.seqdoop.hadoop_bam.SAMFormat;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GaeaAnySAMSortInputFormat extends
		FileInputFormat<LongWritable, SamRecordWritable> {

	public static final String TRUST_EXTS_PROPERTY = "hadoopbam.anysam.trust-exts";
	public static final String SAM_FORMAT_FOR_ALL_PATH = "samformat.allpath";

	private final GaeaBamInputFormat bamIF = new GaeaBamInputFormat();
	private final GaeaSamInputFormat samIF = new GaeaSamInputFormat();

	private final Map<Path, SAMFormat> formatMap;

	private Configuration conf;
	private boolean trustExts;

	public GaeaAnySAMSortInputFormat() {
		this.formatMap = new HashMap<Path, SAMFormat>();
		this.conf = null;
	}

	public GaeaAnySAMSortInputFormat(Configuration conf) {
		this.formatMap = new HashMap<Path, SAMFormat>();
		this.conf = conf;
		this.trustExts = conf.getBoolean(TRUST_EXTS_PROPERTY, true);
	}

	public SAMFormat getFormat(final Path path) {
		SAMFormat fmt = formatMap.get(path);
		if (fmt != null || formatMap.containsKey(path))
			return fmt;
		
		if (this.conf == null)
			throw new IllegalStateException("Don't have a Configuration yet");
		
		if(conf.get(SAM_FORMAT_FOR_ALL_PATH) != null){
			String format = conf.get(SAM_FORMAT_FOR_ALL_PATH);
			if(format.equals("BAM") || format.equals("bam"))
				fmt = SAMFormat.BAM;
			else
				fmt = SAMFormat.SAM;
			formatMap.put(path, fmt);
			return fmt;
		}

		if (trustExts) {
			final SAMFormat f = SAMFormat.inferFromFilePath(path);
			if (f != null) {
				formatMap.put(path, f);
				return f;
			}
		}

		try {
			fmt = SAMFormat.inferFromData(path.getFileSystem(conf).open(path));
		} catch (IOException e) {
		}

		formatMap.put(path, fmt);
		return fmt;
	}


	@Override
	public RecordReader<LongWritable, SamRecordWritable> createRecordReader(
			InputSplit split, TaskAttemptContext ctx)
			throws InterruptedException, IOException {
		final Path path;
		if (split instanceof FileSplit)
			path = ((FileSplit) split).getPath();
		else if (split instanceof FileVirtualSplit)
			path = ((FileVirtualSplit) split).getPath();
		else
			throw new IllegalArgumentException("split '" + split
					+ "' has unknown type: cannot extract path");

		if (this.conf == null)
			this.conf = ctx.getConfiguration();

		final SAMFormat fmt = getFormat(path);
		if (fmt == null)
			throw new IllegalArgumentException(
					"unknown SAM format, cannot create RecordReader: " + path);

		RecordReader<LongWritable, SamRecordWritable> rr;
		switch (fmt) {
		case SAM:
			rr = new GaeaSamSortRecordReader(
					samIF.createRecordReader(split, ctx));
			rr.initialize(split, ctx);
			return rr;
		case BAM:
			rr = new GaeaSamSortRecordReader(
					bamIF.createRecordReader(split, ctx));
			rr.initialize(split, ctx);
			return rr;
		default:
			assert false;
			return null;
		}
	}

	@Override
	public boolean isSplitable(JobContext job, Path path) {
		if (this.conf == null)
			this.conf = job.getConfiguration();

		final SAMFormat fmt = getFormat(path);
		if (fmt == null)
			return super.isSplitable(job, path);

		switch (fmt) {
		case SAM:
			return samIF.isSplitable(job, path);
		case BAM:
			return bamIF.isSplitable(job, path);
		default:
			assert false;
			return false;
		}
	}

	@Override
	public List<InputSplit> getSplits(JobContext job) throws IOException {
		if (this.conf == null)
			this.conf = job.getConfiguration();

		final List<InputSplit> origSplits = super.getSplits(job);

		final List<InputSplit> bamOrigSplits = new ArrayList<InputSplit>(
				origSplits.size()), newSplits = new ArrayList<InputSplit>(
				origSplits.size());

		for (final InputSplit iSplit : origSplits) {
			final FileSplit split = (FileSplit) iSplit;

			if (SAMFormat.BAM.equals(getFormat(split.getPath())))
				bamOrigSplits.add(split);
			else
				newSplits.add(split);
		}
		newSplits.addAll(bamIF.getSplits(bamOrigSplits,
				job.getConfiguration()));
		return newSplits;
	}
}
