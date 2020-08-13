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
package org.bgi.flexlab.gaea.data.mapreduce.input.bam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.tools.mapreduce.bamsort.BamSortOptions;
import org.seqdoop.hadoop_bam.BAMRecordReader;
import org.seqdoop.hadoop_bam.util.MurmurHash3;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;

public class GaeaSamSortRecordReader extends
		RecordReader<LongWritable, SamRecordWritable> {
	public static final String SAMPLENAME_ARRAY_PROP = "bamsort.sample.names";
	private final RecordReader<LongWritable, SamRecordWritable> baseRR;

	private HashMap<String,Integer> sampleID = new HashMap<>();
	private BamSortOptions options = new BamSortOptions();

	public GaeaSamSortRecordReader(RecordReader<LongWritable, SamRecordWritable> rr) {
		baseRR = rr;
	}

	@Override
	public void initialize(InputSplit spl, TaskAttemptContext ctx)
			throws InterruptedException, IOException {
		baseRR.initialize(spl, ctx);

		Configuration conf = ctx.getConfiguration();
		options.getOptionsFromHadoopConf(conf);
		SAMFileHeader header = SamHdfsFileHeader.getHeader(conf);

		List<SAMReadGroupRecord> list = header.getReadGroups();

		for(int i=0;i<list.size();i++)
			sampleID.put(list.get(i).getSample(), i);
	}

	@Override
	public void close() throws IOException {
		if(baseRR != null)
			baseRR.close();
		sampleID.clear();
	}

	@Override
	public float getProgress() throws InterruptedException, IOException {
		return baseRR.getProgress();
	}

	@Override
	public LongWritable getCurrentKey() throws InterruptedException,
			IOException {
		return baseRR.getCurrentKey();
	}

	@Override
	public SamRecordWritable getCurrentValue() throws InterruptedException,
			IOException {
		return baseRR.getCurrentValue();
	}
	
	public long setKey(SAMRecord r) throws InterruptedException, IOException{
		long newKey;
		int ridx = r.getReferenceIndex();
		int start = r.getAlignmentStart();

		String sample = r.getReadGroup().getSample();
		if(!sampleID.containsKey(sample))
			throw new RuntimeException("cantains not sample "+sample+"\t"+sampleID.size());
		int idIndex = sampleID.get(sample);

		if((ridx < 0 || start < 0)){
			int hash = 0;
			byte[] var;
			if ((var = r.getVariableBinaryRepresentation()) != null) {
				// Undecoded BAM record: just hash its raw data.
				hash = (int)MurmurHash3.murmurhash3(var, hash);
			} else {
				// Decoded BAM record or any SAM record: hash a few representative
				// fields together.
				hash = (int)MurmurHash3.murmurhash3(r.getReadName(), hash);
				hash = (int)MurmurHash3.murmurhash3(r.getReadBases(), hash);
				hash = (int)MurmurHash3.murmurhash3(r.getBaseQualities(), hash);
				hash = (int)MurmurHash3.murmurhash3(r.getCigarString(), hash);
			}
			hash = Math.abs(hash);
			newKey = ((long)idIndex << 48) | ((long)65535 << 32) | (long)hash;
		}else{
			newKey = ((long)idIndex << 48) | (((long)ridx) << 32) | (long)start;
		}
		getCurrentKey().set(newKey);

		return  newKey;
	}

	@Override
	public boolean nextKeyValue() throws InterruptedException, IOException {
		if (!baseRR.nextKeyValue()){
			return false;
		}

		final SAMRecord r = getCurrentValue().get();
		setKey(r);
		return true;
	}
}
