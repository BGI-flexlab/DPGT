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

import htsjdk.samtools.SAMRecord;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.WritableComparable;
import org.bgi.flexlab.gaea.data.exception.OutOfBoundException;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.HashMap;

public class WindowsBasedWritable implements WritableComparable<WindowsBasedWritable> {
	private LongWritable windowsInfo = new LongWritable();
	private IntWritable position = new IntWritable();

	private final static int SAMPLE_BITS = 22;
	private final static int WINDOW_NUMBER_BITS = 32;
	private final static int CHROMOSOME_BITS = Long.BYTES * Byte.SIZE - SAMPLE_BITS - WINDOW_NUMBER_BITS;

	private final static int MAX_SAMPLE_ID = (int) (Math.pow(2, SAMPLE_BITS));
	private final static int CHROMOSOME_BIT_INDEX = SAMPLE_BITS + WINDOW_NUMBER_BITS;

	private final static int CHROMOSOME_BITS_MASK = (int) (Math.pow(2, CHROMOSOME_BITS) - 1);
	private final static int SAMPLE_BITS_MASK = (int) (Math.pow(2, SAMPLE_BITS) - 1);
	private final static int WINDOW_NUMBER_MASK = (int) (Math.pow(2, WINDOW_NUMBER_BITS) - 1);

	//added by gc
	//是为了对位置进行全局排序，进而让group有序的分配给reduce，而不是穿插的顺序
	public static HashMap<Integer,Long> accumulate=new HashMap<>();
	public static Long realStart=0L;
	public static Long realEnd=0L;
	public static Integer partitionsNum=1;
	public WindowsBasedWritable(HashMap<Integer,Long> chrAccLength,Long start,Long end,Integer parNum){
		accumulate=chrAccLength;
		realStart=start;
		realEnd=end;
		partitionsNum=parNum;
	}
	public WindowsBasedWritable(){
	}
	public void set(long sample, long chromosome, long winNum, int pos) {
		if (sample >= MAX_SAMPLE_ID)
			throw new OutOfBoundException(String.format("sample size %d is more than 4194304", (int) sample));

		if (winNum >= Integer.MAX_VALUE)
			throw new OutOfBoundException(
					String.format("window number %d is more than %d", (int) winNum, Integer.MAX_VALUE));

		if (chromosome >= CHROMOSOME_BITS_MASK)
			throw new OutOfBoundException(
					String.format("chromosome size %d is more than %d", (int) chromosome, CHROMOSOME_BITS_MASK));

		long key = 0;
		key = (chromosome << CHROMOSOME_BIT_INDEX) | (winNum << SAMPLE_BITS) | sample;
		windowsInfo.set(key);
		position.set(pos);
	}

	public void set(int chromosome, int winNum, int pos) {
		set(0, chromosome, winNum, pos);
	}

	public String toString() {
		return windowsInfo.toString() + "\t" + position.get();
	}

	public int getChromosomeIndex() {
		long key = windowsInfo.get();
		int index = (int) ((key >> CHROMOSOME_BIT_INDEX) & CHROMOSOME_BITS_MASK);
		if (index == CHROMOSOME_BITS_MASK)
			return SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
		return index;
	}

	public long getWindows() {
		return windowsInfo.get();
	}

	public LongWritable getWindowsInformation() {
		return windowsInfo;
	}

	public int getWindowsNumber() {
		long key = windowsInfo.get();
		return (int) ((key >> SAMPLE_BITS) & WINDOW_NUMBER_MASK);
	}

	public int getSampleID() {
		long key = windowsInfo.get();
		return (int) (key & SAMPLE_BITS_MASK);
	}

	public IntWritable getPosition() {
		return position;
	}

	@Override
	public void readFields(DataInput in) throws IOException {
		windowsInfo.readFields(in);
		position.readFields(in);
	}

	public void write(DataOutput out) throws IOException {
		windowsInfo.write(out);
		position.write(out);
	}

	@Override
	public int hashCode() {
		return ((windowsInfo.hashCode() * 163 + position.hashCode()) & 0xffffffff);
	}

	//修改成全局排序的

	public int partition() {
		/*原代码
		int hashcode = (getChromosomeIndex() + 1);
		hashcode += (getWindowsNumber()+1);
		hashcode += (getSampleID() + 1);
		return (int)(hashcode & 0xffffffff);
		*/

		//修改后代码
		if(accumulate.size()>0) {
			//全局排序方式
			Integer chrIndex = getChromosomeIndex();
			Long curPos = position.get() + accumulate.get(chrIndex);
			Long band = (long) ((realEnd - realStart) * 1.1 / partitionsNum);
			if (band == 0) {
				System.out.println("input data is empty or un-recognized");
				System.exit(-1);
			}
			return (int) ((curPos - realStart) / band);
		}else{
			//老方式，穿插的排序
			int hashcode = (getChromosomeIndex() + 1);
			hashcode += (getWindowsNumber()+1);
			hashcode += (getSampleID() + 1);
			return (int)(hashcode & 0xffffffff);
		}
	}

	@Override
	public boolean equals(Object other) {
		if (other instanceof WindowsBasedWritable) {
			WindowsBasedWritable tmp = (WindowsBasedWritable) other;
			return windowsInfo.equals(tmp.getWindowsInformation())
					&& position.get() == (tmp.getPosition().get());
		}
		return false;
	}

	@Override
	public int compareTo(WindowsBasedWritable tp) {
		long cmp = windowsInfo.get() - tp.getWindows();
		if (cmp != 0) {
			if (cmp > 0)
				return 1;
			return -1;
		}
		return position.compareTo(tp.position);
	}
}
