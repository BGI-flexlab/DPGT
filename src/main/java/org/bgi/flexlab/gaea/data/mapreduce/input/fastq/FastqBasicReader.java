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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.CompressionCodecFactory;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.apache.hadoop.util.LineReader;

import java.io.IOException;

public abstract class FastqBasicReader {
	protected static final Log LOG = LogFactory.getLog(FastqBasicReader.class
			.getName());

	protected CompressionCodecFactory compressionCodecs = null;
	protected long start;
	protected long pos;
	protected long end;
	protected LineReader in;
	protected int maxLineLength;
	protected byte[] recordDelimiterBytes;
	protected String firstLine = "";
	protected String sampleID;
	protected String secondLine = "";

	public void getFirstFastqLine() throws IOException {
		Text tmpline = new Text();
		int size;
		while ((size = in.readLine(tmpline, maxLineLength, Math.max(
				(int) Math.min(Integer.MAX_VALUE, end - pos), maxLineLength))) != 0) {
			start += size;
			if (tmpline.toString().startsWith("@")) {
				firstLine = tmpline.toString();
				if ((size = in.readLine(tmpline, maxLineLength, Math.max(
						(int) Math.min(Integer.MAX_VALUE, end - pos),
						maxLineLength))) != 0) {
					start += size;
					if (tmpline.toString().startsWith("@")) {
						firstLine = tmpline.toString();
						if ((size = in.readLine(tmpline,maxLineLength,
								Math.max((int) Math.min(Integer.MAX_VALUE, end
												- pos), maxLineLength))) != 0) {
							start += size;
							secondLine = tmpline.toString();
						}
					} else {
						secondLine = tmpline.toString();
					}
				}
				break;
			}
		}
	}

	public void setStart(long start) {
		this.start = start;
	}

	public void setEnd(long end) {
		this.end = end;
	}

	public FastqBasicReader(Configuration job, FileSplit split,
			byte[] recordDelimiter) throws IOException {
		this.maxLineLength = job.getInt("mapred.linerecordreader.maxlength",
				Integer.MAX_VALUE);
		compressionCodecs = new CompressionCodecFactory(job);
		final CompressionCodec codec = compressionCodecs.getCodec(split
				.getPath());

		String multiSampleList = job.get("multiSampleList");
		if (multiSampleList != null && multiSampleList != "") {
			FastqMultipleSample samplelist;
			samplelist = new FastqMultipleSample(multiSampleList, false);
			FastqSample slist = samplelist.getID(split.getPath().toString());
			if (slist != null) {
				sampleID = String.valueOf(slist.getId());
			} else {
				sampleID = "+";
			}
		}

		start = split.getStart();
		end = split.getStart() + split.getLength();

		// open the file and seek to the start of the split
		FileSystem fs = split.getPath().getFileSystem(job);
		FSDataInputStream fileIn = fs.open(split.getPath());
		boolean skipFirstLine = false;
		if (codec != null) {
			if (null == this.recordDelimiterBytes) {
				in = new LineReader(codec.createInputStream(fileIn), job);
			} else {
				in = new LineReader(codec.createInputStream(fileIn), job,
						this.recordDelimiterBytes);
			}
			end = Long.MAX_VALUE;
		} else {
			if (start != 0) {
				skipFirstLine = true;
				--start;
				fileIn.seek(start);
			}
			if (null == this.recordDelimiterBytes) {
				in = new LineReader(fileIn, job);
			} else {
				in = new LineReader(fileIn, job, this.recordDelimiterBytes);
			}
		}

		if (skipFirstLine) { // skip first line and re-establish "start".
			start += in.readLine(new Text(), 0,
					(int) Math.min((long) Integer.MAX_VALUE, end - start));
		}
		getFirstFastqLine();
		this.pos = start;
	}

	/**
	 * Get the progress within the split
	 */
	public float getProgress() {
		if (start == end) {
			return 0.0f;
		} else {
			return Math.min(1.0f, (pos - start) / (float) (end - start));
		}
	}

	public synchronized void close() throws IOException {
		if (in != null) {
			in.close();
		}
	}

	public synchronized long getPos() throws IOException {
		return pos;
	}

	public abstract boolean next(Text key, Text value) throws IOException;
}
