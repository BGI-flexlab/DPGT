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
package org.bgi.flexlab.gaea.data.mapreduce.writable;

import org.apache.hadoop.io.WritableComparable;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasicCodec;
import org.seqdoop.hadoop_bam.util.DataInputWrapper;
import org.seqdoop.hadoop_bam.util.DataOutputWrapper;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public class AlignmentBasicWritable implements WritableComparable<AlignmentBasicWritable> {
	private static final AlignmentsBasicCodec codec = new AlignmentsBasicCodec();

	private AlignmentsBasic alignment;

	public AlignmentsBasic getAlignment() {
		return alignment;
	}

	public void setAlignment(AlignmentsBasic alignment) {
		this.alignment = alignment;
	}

	@Override
	public void readFields(DataInput dataInput) throws IOException {
		codec.setInputStream(new DataInputWrapper(dataInput));
		alignment = codec.decode();
	}

	@Override
	public void write(DataOutput dataOutput) throws IOException {
		//final AlignmentsBasicCodec codec = new AlignmentsBasicCodec();
		codec.setOutputStream(new DataOutputWrapper(dataOutput));
		codec.encode(alignment);
	}

	@Override
	public int compareTo(AlignmentBasicWritable o) {
		AlignmentsBasic alignment2 = o.getAlignment();
		if(alignment.getChrNameIndex() != alignment2.getChrNameIndex()) {
			if(alignment.getChrNameIndex() > alignment2.getChrNameIndex()) {
				return 1;
			} else {
				return -1;
			}
		} else {
			if(alignment.getPosition() != alignment2.getPosition()) {
				if(alignment.getPosition() > alignment2.getPosition()) {
					return 1;
				} else {
					return -1;
				}
			}
		}
		return 0;
	}
	
}