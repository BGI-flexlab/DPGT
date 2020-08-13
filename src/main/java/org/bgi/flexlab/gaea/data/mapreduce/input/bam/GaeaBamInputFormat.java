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

import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.seqdoop.hadoop_bam.FileVirtualSplit;
import org.seqdoop.hadoop_bam.SplittingBAMIndex;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class GaeaBamInputFormat extends
		FileInputFormat<LongWritable, SamRecordWritable> {
	public static boolean DEBUG_BAM_SPLITTER = false; 

	private Path getIdxPath(Path path) {
		return path.suffix(".splitting-bai");
	}

	public RecordReader<LongWritable, SamRecordWritable> createRecordReader(
			InputSplit split, TaskAttemptContext ctx)
			throws InterruptedException, IOException {
		RecordReader<LongWritable, SamRecordWritable> rr = new GaeaBamRecordReader();
		Configuration conf = ctx.getConfiguration();
		DEBUG_BAM_SPLITTER = conf.getBoolean("debug.bam.splitter", false);
		rr.initialize(split, ctx);
		return rr;
	}

	public List<InputSplit> getSplits(JobContext job) throws IOException {
		return getSplits(super.getSplits(job), job.getConfiguration());
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public List<InputSplit> getSplits(List<InputSplit> splits, Configuration cfg)
			throws IOException {
		Collections.sort(splits, new Comparator() {
			@SuppressWarnings("unused")
			public int compare(InputSplit a, InputSplit b) {
				FileSplit fa = (FileSplit) a;
				FileSplit fb = (FileSplit) b;
				return fa.getPath().compareTo(fb.getPath());
			}

			@Override
			public int compare(Object a, Object b) {
				FileSplit fa = (FileSplit) a;
				FileSplit fb = (FileSplit) b;
				return fa.getPath().compareTo(fb.getPath());
			}
		});
		List<InputSplit> newSplits = new ArrayList<InputSplit>(splits.size());

		for (int i = 0; i < splits.size();) {
			try {
				i = addIndexedSplits(splits, i, newSplits, cfg);
			} catch (IOException e) {
				i = addProbabilisticSplits(splits, i, newSplits, cfg);
			}
		}
		return newSplits;
	}

	private int addIndexedSplits(List<InputSplit> splits, int i,
			List<InputSplit> newSplits, Configuration cfg) throws IOException {
		Path file = ((FileSplit) splits.get(i)).getPath();

		SplittingBAMIndex idx = new SplittingBAMIndex(file.getFileSystem(cfg)
				.open(getIdxPath(file)));

		int splitsEnd = splits.size();
		for (int j = i; j < splitsEnd; j++) {
			if (!file.equals(((FileSplit) splits.get(j)).getPath()))
				splitsEnd = j;
		}
		for (int j = i; j < splitsEnd; j++) {
			FileSplit fileSplit = (FileSplit) splits.get(j);

			long start = fileSplit.getStart();
			long end = start + fileSplit.getLength();

			Long blockStart = idx.nextAlignment(start);

			Long blockEnd = Long.valueOf(j == splitsEnd - 1 ? idx
					.prevAlignment(end).longValue() | 0xFFFF : idx
					.nextAlignment(end).longValue());

			if (blockStart == null) {
				throw new RuntimeException(
						"Internal error or invalid index: no block start for "
								+ start);
			}
			if (blockEnd == null) {
				throw new RuntimeException(
						"Internal error or invalid index: no block end for "
								+ end);
			}
			newSplits.add(new FileVirtualSplit(file, blockStart.longValue(),
					blockEnd.longValue(), fileSplit.getLocations()));
		}
		return splitsEnd;
	}

	private int addProbabilisticSplits(List<InputSplit> splits, int i,
			List<InputSplit> newSplits, Configuration cfg) throws IOException {
		Path path = ((FileSplit) splits.get(i)).getPath();
		SeekableStream sin = WrapSeekable.openPath(path.getFileSystem(cfg),
				path);

		GaeaBamSplitGuesser guesser = new GaeaBamSplitGuesser(sin,cfg);

		FileVirtualSplit previousSplit = null;

		for (; i < splits.size(); i++) {
			FileSplit fspl = (FileSplit) splits.get(i);
			if (!fspl.getPath().equals(path)) {
				break;
			}
			long beg = fspl.getStart();
			long end = beg + fspl.getLength();

			long alignedBeg = guesser.guessNextBAMRecordStart(beg, end);

			long alignedEnd = end << 16 | 0xFFFF;

			if (alignedBeg == end) {
				if (previousSplit == null) {
					System.err
							.println("'"
									+ path
									+ "': "
									+ "no reads in first split: bad BAM file or tiny split size?");
				} else {
					previousSplit.setEndVirtualOffset(alignedEnd);
				}
			} else
				newSplits.add(previousSplit = new FileVirtualSplit(path,
						alignedBeg, alignedEnd, fspl.getLocations()));

		}

		sin.close();
		return i;
	}

	public boolean isSplitable(JobContext job, Path path) {
		return true;
	}

}
