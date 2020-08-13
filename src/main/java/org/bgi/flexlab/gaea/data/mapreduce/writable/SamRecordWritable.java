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
package org.bgi.flexlab.gaea.data.mapreduce.writable;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMRecord;
import org.apache.hadoop.io.Writable;
import org.seqdoop.hadoop_bam.LazyBAMRecordFactory;
import org.seqdoop.hadoop_bam.util.DataInputWrapper;
import org.seqdoop.hadoop_bam.util.DataOutputWrapper;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public class SamRecordWritable implements Writable{
	private static final BamRecordCodec lazyCodec =
			new BamRecordCodec(null, new LazyBAMRecordFactory());

		private SAMRecord record;

		public SAMRecord get()            { return record; }
		public void      set(SAMRecord r) { record = r; }

		@Override 
		public void write(DataOutput out) throws IOException {
			// In theory, it shouldn't matter whether we give a header to
			// BAMRecordCodec or not, since the representation of an alignment in BAM
			// doesn't depend on the header data at all. Only its interpretation
			// does, and a simple read/write codec shouldn't really have anything to
			// say about that. (But in practice, it already does matter for decode(),
			// which is why LazyBAMRecordFactory exists.)
			final BAMRecordCodec codec = new BAMRecordCodec(record.getHeader());
			codec.setOutputStream(new DataOutputWrapper(out));
			codec.encode(record);
		}
		
		@Override 
		public void readFields(DataInput in) throws IOException {
			lazyCodec.setInputStream(new DataInputWrapper(in));
			record = lazyCodec.decode();
		}

		@Override
		public String toString() {
			return record.getSAMString().trim(); // remove trailing newline
		}
}
