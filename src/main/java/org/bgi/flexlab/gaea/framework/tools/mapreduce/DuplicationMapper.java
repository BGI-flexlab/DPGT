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
package org.bgi.flexlab.gaea.framework.tools.mapreduce;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.Writable;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.CreateDuplicationKey;
import org.bgi.flexlab.gaea.data.mapreduce.writable.DuplicationKeyWritable;
import org.seqdoop.hadoop_bam.SAMRecordWritable;

public class DuplicationMapper extends PairEndAggregatorMapper{
	
	protected SAMFileHeader header = null;
	protected CreateDuplicationKey dupKey = null;
	protected DuplicationKeyWritable dupKeyWritable = null;
	
	@Override
	public void setup(Context context) {
		Configuration conf = context.getConfiguration();
		header = SamHdfsFileHeader.getHeader(conf);
		dupKey = new CreateDuplicationKey(header);
		dupKeyWritable = new DuplicationKeyWritable();
	}
	
	protected Writable getKey(Writable keyin,Writable valuein){
		if(valuein instanceof SAMRecordWritable){
			SAMRecord sam = ((SAMRecordWritable)valuein).get();
			dupKey.getKey(sam,dupKeyWritable);
		}
		return dupKeyWritable;
	}
}
