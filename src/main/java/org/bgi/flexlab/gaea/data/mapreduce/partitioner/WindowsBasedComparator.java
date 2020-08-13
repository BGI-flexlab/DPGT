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

import org.apache.hadoop.io.DataInputBuffer;
import org.apache.hadoop.io.RawComparator;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;

import java.io.IOException;

public class WindowsBasedComparator implements RawComparator<WindowsBasedWritable> {

	@Override
	public int compare(WindowsBasedWritable o1, WindowsBasedWritable o2) {
		if(o1.getWindows() == o2.getWindows())
			return 0;
		if(o1.getWindows() > o2.getWindows())
			return 1;
		return -1;
	}

	@Override
	public int compare(byte[] b1, int s1, int l1, byte[] b2, int s2, int l2) {
		WindowsBasedWritable key1 = new WindowsBasedWritable();
		WindowsBasedWritable key2 = new WindowsBasedWritable();
		DataInputBuffer buffer = new DataInputBuffer();
		try {
			buffer.reset(b1, s1, l1);
			key1.readFields(buffer);
			buffer.reset(b2, s2, l2);
			key2.readFields(buffer);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return compare(key1, key2);
	}
}
