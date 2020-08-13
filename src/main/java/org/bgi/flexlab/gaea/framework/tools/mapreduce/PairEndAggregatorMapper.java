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

import htsjdk.samtools.SAMRecord;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.hadoop.mapreduce.Mapper;
import org.seqdoop.hadoop_bam.SAMRecordWritable;

import java.io.IOException;

public class PairEndAggregatorMapper extends
		Mapper<Writable, Writable, Writable, Writable> {
	protected Writable keyout ;
	protected Text textKey = new Text();
	protected Writable valueout;
	
	protected Writable getKey(Writable keyin,Writable valuein){
		if(keyin instanceof Text){
			textKey.set((Text)keyin);
		}else if(keyin instanceof LongWritable){
			SAMRecord sam = ((SAMRecordWritable)valuein).get();
			textKey.set(sam.getReadName());
		}
		return textKey;
	}
	
	protected Writable getValue(Writable valuein){
		return valuein;
	}

	protected void map(Writable key, Writable value, Context context)
			throws IOException, InterruptedException {
		keyout = getKey(key,value);
		valueout = getValue(value);
		context.write(keyout, valueout);
	}
}
