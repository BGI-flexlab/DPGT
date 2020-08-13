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
package org.bgi.flexlab.gaea.data.mapreduce.output.cram;

import htsjdk.samtools.CRAMContainerStreamWriter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTagUtil;
import htsjdk.samtools.cram.ref.ReferenceSource;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.seqdoop.hadoop_bam.SAMRecordWritable;

import java.io.IOException;
import java.io.OutputStream;
import java.net.URI;
import java.nio.file.Paths;

public class GaeaKeyIgnoringCramRecordWriter<K> extends
		RecordWriter<K, SAMRecordWritable> {
	public final static String OUTPUTFORMAT_REFERENCE = "cram.outputformat.reference";
	private static final String HADOOP_BAM_PART_ID = "Hadoop-BAM-Part";
	private OutputStream origOutput;
	private CRAMContainerStreamWriter cramContainerStream = null;
	private ReferenceSource refSource = null;
	private boolean writeHeader = true;
	private boolean rename = false;
	private String sample = null;
	private Path outputPath = null;

	public GaeaKeyIgnoringCramRecordWriter(Path path, Boolean writerHeader,
			TaskAttemptContext ctx) throws IOException {
		this.outputPath = path;
		init(path, writerHeader, ctx);
	}

	private void init(final Path output, final boolean writeHeader,
			final TaskAttemptContext ctx) throws IOException {
		init(output.getFileSystem(ctx.getConfiguration()).create(output),
				writeHeader, ctx);
	}

	private void init(final OutputStream output, final boolean writeHeader,
			final TaskAttemptContext ctx) throws IOException {
		origOutput = output;
		this.writeHeader = writeHeader;

		final URI referenceURI = URI.create(ctx.getConfiguration().get(
				OUTPUTFORMAT_REFERENCE));

		rename = ctx.getConfiguration().getBoolean("rename.file", false);

		refSource = new ReferenceSource(Paths.get(referenceURI));
	}

	@Override
	public void close(TaskAttemptContext ctx) throws IOException {
		cramContainerStream.finish(true);
		origOutput.close();

		if (rename) {
			final FileSystem srcFS = outputPath.getFileSystem(ctx
					.getConfiguration());
			if (this.sample != null) {
				Path newName = new Path(outputPath.getParent() + "/" + sample
						+ ".sorted.cram");
				srcFS.rename(outputPath, newName);
			}
		}
	}

	protected void writeAlignment(final SAMRecord rec) {
		if (null == cramContainerStream) {
			String rgId = (String) rec
					.getAttribute(SAMTagUtil.getSingleton().RG);
			final SAMFileHeader header = rec.getHeader();
			if (header == null) {
				throw new RuntimeException(
						"Cannot write record to CRAM: null header in SAM record");
			}
			sample = header.getReadGroup(rgId).getSample();
			if (writeHeader) {
				this.writeHeader(header);
			}
			cramContainerStream = new CRAMContainerStreamWriter(origOutput,
					null, refSource, header, HADOOP_BAM_PART_ID);
		}
		cramContainerStream.writeAlignment(rec);
	}

	private void writeHeader(final SAMFileHeader header) {
		cramContainerStream.writeHeader(header);
	}

	@Override
	public void write(K arg0, SAMRecordWritable writable) throws IOException,
			InterruptedException {
		writeAlignment(writable.get());
	}
}
