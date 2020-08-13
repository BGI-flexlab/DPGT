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

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public class PairWritable implements WritableComparable<PairWritable> {
    private String first;
    private String second;

    public PairWritable() {
    }

    public PairWritable(String first, String second) {
        this.set(first, second);
    }

    public void set(String first, String second) {
        this.first = first;
        this.second = second;
    }

    /**
     * 反序列化
     */
    @Override
    public void readFields(DataInput arg0) throws IOException {
        this.first = arg0.readUTF();
        this.second = arg0.readUTF();
    }

    @Override
    public void write(DataOutput arg0) throws IOException {
        arg0.writeUTF(first);
        arg0.writeUTF(second);
    }

    @Override
    public int compareTo(PairWritable o) {
        int comp = this.first.compareTo(o.first);
        return comp != 0 ? comp : this.second.compareTo(o.getSecond());
    }

    public String getSecond() {
        return second;
    }

    public void setSecond(String second) {
        this.second = second;
    }

    public String getFirst() {
        return first;
    }

    public void setFirst(String first) {
        this.first = first;
    }
}