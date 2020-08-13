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
package org.bgi.flexlab.gaea.data.mapreduce.partitioner;


import org.apache.hadoop.mapreduce.Partitioner;
import org.bgi.flexlab.gaea.data.mapreduce.writable.PairWritable;

public class FirstPartitioner<T> extends Partitioner<PairWritable, T> {

    @Override
    public int getPartition(PairWritable key, T value, int numPartitions) {
        /*
         * 默认的实现 (key.hashCode() & Integer.MAX_VALUE) % numPartitions
         * 让key中first字段作为分区依据
         */
        return (key.getFirst().hashCode() & Integer.MAX_VALUE) % numPartitions;
    }
}
