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
package org.bgi.flexlab.gaea.data.mapreduce.output.bam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;

import java.io.IOException;
import java.io.OutputStream;
import java.io.StringWriter;
import java.io.Writer;

public class GaeaKeyIgnoringBamRecordWriter<K> extends
		RecordWriter<K, SamRecordWritable> {

	protected SAMFileHeader header;
	private BinaryCodec binaryCodec;
	private BAMRecordCodec bamRecordCodec;
	private Path outputPath;
	private OutputStream outputStream;
	private boolean writeHeader;

	public GaeaKeyIgnoringBamRecordWriter(Path p, Boolean w,
			TaskAttemptContext ctx) throws IOException {
		this.outputPath = p;
		this.outputStream = outputPath.getFileSystem(ctx.getConfiguration()).create(
				outputPath);
		this.writeHeader = w;
	}

	public GaeaKeyIgnoringBamRecordWriter(Path p, SAMFileHeader header,Boolean w,
										  TaskAttemptContext ctx) throws IOException {
		this.outputPath = p;
		this.outputStream = outputPath.getFileSystem(ctx.getConfiguration()).create(
				outputPath);
		this.writeHeader = w;
		this.header = header;
		initialize(header);
	}

	public GaeaKeyIgnoringBamRecordWriter(OutputStream os, Boolean w,
			TaskAttemptContext ctx) {
		this.outputStream = os;
		this.writeHeader = w;
	}

	private void initialize(SAMFileHeader header) {
		OutputStream compressedOut = null;
		if (outputStream != null)
			compressedOut = new BlockCompressedOutputStream(outputStream, null);
		else
			compressedOut = new BlockCompressedOutputStream(outputPath.toString());

		binaryCodec = new BinaryCodec(compressedOut);
		bamRecordCodec = new BAMRecordCodec(header);
		bamRecordCodec.setOutputStream(compressedOut);

		if (writeHeader) {
			writeHeader(header);
		}
	}

	private void writeHeader(final SAMFileHeader header) {
		System.err.println(outputPath.toString() + "\t" + outputStream.toString());
		binaryCodec.writeBytes("BAM\001".getBytes());

		final Writer sw = new StringWriter();
		new SAMTextHeaderCodec().encode(sw, header);

		binaryCodec.writeString(sw.toString(), true, false);

		final SAMSequenceDictionary dict = header.getSequenceDictionary();

		binaryCodec.writeInt(dict.size());
		for (final SAMSequenceRecord rec : dict.getSequences()) {
			binaryCodec.writeString(rec.getSequenceName(), true, true);
			binaryCodec.writeInt(rec.getSequenceLength());
		}
	}

	@Override
	public void close(TaskAttemptContext context) throws IOException,
			InterruptedException {
		if (binaryCodec != null) {
			binaryCodec.close();
		}
	}

	@Override
	public void write(K ignored, SamRecordWritable samRecordWritable) throws IOException,
			InterruptedException {
		SAMRecord sam = samRecordWritable.get();
		if (binaryCodec == null || bamRecordCodec == null) {
			initialize(sam.getHeader());
		}
		
//		/*check?*/
//		if(sam.getReadUnmappedFlag()){
//			sam.setAlignmentStart(0);
//		}
		bamRecordCodec.encode(sam);
	}
}
