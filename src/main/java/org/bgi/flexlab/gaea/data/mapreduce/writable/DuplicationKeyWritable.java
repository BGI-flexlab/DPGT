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

import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.hadoop.io.BooleanWritable;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.WritableComparable;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public class DuplicationKeyWritable implements
		WritableComparable<DuplicationKeyWritable> {
	private Text LB;
	private IntWritable chrIndex;
	private IntWritable position;
	private BooleanWritable forward;

	public DuplicationKeyWritable() {
		LB = new Text();
		chrIndex = new IntWritable();
		position = new IntWritable();
		forward = new BooleanWritable();
	}

	public DuplicationKeyWritable(Text LB, IntWritable chrIndex, IntWritable position,
			BooleanWritable forward) {
		this.LB = LB;
		this.chrIndex = chrIndex;
		this.position = position;
		this.forward = forward;
	}

	public void set(String LB, int chrIndex, int position, boolean forward) {
		this.LB.set(LB);
		this.chrIndex.set(chrIndex);
		this.position.set(position);
		this.forward.set(forward);
	}

	@Override
	public void readFields(DataInput in) throws IOException {
		LB.readFields(in);
		chrIndex.readFields(in);
		position.readFields(in);
		forward.readFields(in);
	}

	@Override
	public void write(DataOutput out) throws IOException {
		LB.write(out);
		chrIndex.write(out);
		position.write(out);
		forward.write(out);
	}

	@Override
	public int compareTo(DuplicationKeyWritable key) {
		int v1 = LB.compareTo(key.LB);
		int v2 = chrIndex.compareTo(key.chrIndex);
		int v3 = position.compareTo(key.position);
		if (v1 != 0)
			return v1;
		if (v2 != 0)
			return v2;
		if (v3 != 0)
			return v3;
		if (forward.compareTo(key.forward) != 0) {
			if (forward.get())
				return 1;
			else
				return -1;
		}
		return 0;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder().append(LB).append(chrIndex)
				.append(position).append(forward).toHashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof DuplicationKeyWritable)) {
			return false;
		} else if (obj == this)
			return true;

		DuplicationKeyWritable key = (DuplicationKeyWritable) obj;
		if (this.compareTo(key) == 0)
			return true;

		return false;
	}

	public String getLB() {
		return LB.toString();
	}

	public int getChrIndex() {
		return chrIndex.get();
	}

	public int getPosition() {
		return position.get();
	}

	public boolean isForward() {
		return forward.get();
	}
}
