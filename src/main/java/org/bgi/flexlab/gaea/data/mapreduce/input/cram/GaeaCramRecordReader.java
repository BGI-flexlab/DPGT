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

import htsjdk.samtools.CRAMIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.seqdoop.hadoop_bam.util.MurmurHash3;
import org.seqdoop.hadoop_bam.util.SAMHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import java.io.File;
import java.io.IOException;

public class GaeaCramRecordReader extends
		RecordReader<LongWritable, SamRecordWritable> {
	public final static String INPUTFORMAT_REFERENCE = "inputformat.reference";
	public final static String CRAM_FILE_SPLITABLE = "cram.file.splitable";

	protected final LongWritable key = new LongWritable();
	protected final SamRecordWritable record = new SamRecordWritable();
	private boolean isInitialized = false;
	protected SeekableStream seekableStream;
	protected long start;
	protected long length;
	private CRAMIterator cramIterator;
	protected SAMFileHeader samFileHeader = null;

	@Override
	public void initialize(InputSplit split, TaskAttemptContext context)
			throws IOException {
		if (isInitialized) {
			close();
		}
		isInitialized = true;

		final Configuration conf = context.getConfiguration();
		final FileSplit fileSplit = (FileSplit) split;
		final Path file = fileSplit.getPath();

		String refSourcePath = conf.get(INPUTFORMAT_REFERENCE);

		ReferenceSource refSource = new ReferenceSource(new File(refSourcePath));
		
		seekableStream = WrapSeekable.openPath(conf, file);
		start = getStart(fileSplit, conf);
		if (start == 0) {
			samFileHeader = CramIO.readCramHeader(seekableStream).getSamFileHeader();
			start = seekableStream.position();
			seekableStream.seek(0);
		}

		length = getLength(fileSplit, conf, seekableStream.length());
		long end = start + length;
		if (end > seekableStream.length())
			end = seekableStream.length();

		long[] boundaries = new long[] { start << 16, (end - 1) << 16 };
		cramIterator = new CRAMIterator(seekableStream, refSource, boundaries,ValidationStringency.DEFAULT_STRINGENCY);
		ValidationStringency stringency = SAMHeaderReader
				.getValidationStringency(conf);
		if (stringency != null) {
			cramIterator.setValidationStringency(stringency);
		}
	}

	private long getStart(FileSplit fileSplit, Configuration conf) {
		boolean isSplitable = conf.getBoolean(CRAM_FILE_SPLITABLE, false);
		if (isSplitable)
			return fileSplit.getStart();
		return 0;
	}

	private long getLength(FileSplit fileSplit, Configuration conf,
			long fileLength) {
		boolean isSplitable = conf.getBoolean(CRAM_FILE_SPLITABLE, false);
		if (isSplitable)
			return fileSplit.getLength();
		return fileLength;
	}

	@Override
	public boolean nextKeyValue() {
		if (!cramIterator.hasNext()) {
			return false;
		}
		SAMRecord r = cramIterator.next();
		key.set(getKey(r));
		r.setHeader(samFileHeader);
		record.set(r);
		return true;
	}

	@Override
	public LongWritable getCurrentKey() {
		return key;
	}

	@Override
	public SamRecordWritable getCurrentValue() {
		return record;
	}

	@Override
	public float getProgress() throws IOException {
		return (float) (seekableStream.position() - start) / length;
	}

	@Override
	public void close() {
		cramIterator.close();
	}

	/**
	 * Note: this is the only getKey function that handles unmapped reads
	 * specially!
	 */
	private long getKey(final SAMRecord rec) {
		final int refIdx = rec.getReferenceIndex();
		final int start = rec.getAlignmentStart();

		if (!(rec.getReadUnmappedFlag() || refIdx < 0 || start < 0))
			return getKey(refIdx, start);

		int hash = 0;
		byte[] var;
		if ((var = rec.getVariableBinaryRepresentation()) != null) {
			hash = (int) MurmurHash3.murmurhash3(var, hash);
		} else {
			hash = (int) MurmurHash3.murmurhash3(rec.getReadName(), hash);
			hash = (int) MurmurHash3.murmurhash3(rec.getReadBases(), hash);
			hash = (int) MurmurHash3.murmurhash3(rec.getBaseQualities(), hash);
			hash = (int) MurmurHash3.murmurhash3(rec.getCigarString(), hash);
		}
		hash = Math.abs(hash);
		return getKey0(Integer.MAX_VALUE, hash);
	}

	/**
	 * @param alignmentStart
	 *            1-based leftmost coordinate.
	 */
	private long getKey(int refIdx, int alignmentStart) {
		return getKey0(refIdx, alignmentStart - 1);
	}

	/**
	 * @param alignmentStart0
	 *            0-based leftmost coordinate.
	 */
	private long getKey0(int refIdx, int alignmentStart0) {
		return (long) refIdx << 32 | alignmentStart0;
	}
}
