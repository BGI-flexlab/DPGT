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

import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.WritableComparable;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public class WindowsBasedBasicWritable implements WritableComparable<WindowsBasedBasicWritable> {
	protected Text windowsInfo = new Text();

	public void set(String sample, String chromosome, int winNum, int pos) {
		set(sample + ":" + chromosome + ":" + winNum, pos);
	}

	public void set(String chromosome, int winNum, int pos) {
		set(chromosome + ":" + winNum, pos);
	}

	public void set(String winInfo, int pos) {
		this.windowsInfo.set(winInfo);
	}

	public String toString() {
		return windowsInfo.toString();
	}

	public String getChromosomeName() {
		String[] win = windowsInfo.toString().split(":");
		return win[win.length - 2];
	}

	public Text getWindows() {
		return windowsInfo;
	}

	public String getWindowsInformation() {
		return windowsInfo.toString();
	}

	public int getWindowsNumber() {
		String[] win = windowsInfo.toString().split(":");
		return Integer.parseInt(win[win.length - 1]);
	}

	@Override
	public void readFields(DataInput in) throws IOException {
		windowsInfo.readFields(in);
	}

	public void write(DataOutput out) throws IOException {
		windowsInfo.write(out);
	}

	@Override
	public int hashCode() {
		return windowsInfo.hashCode() * 163;
	}

	@Override
	public boolean equals(Object other) {
		if (other instanceof WindowsBasedBasicWritable) {
			WindowsBasedBasicWritable tmp = (WindowsBasedBasicWritable) other;
			return windowsInfo.toString().equals(tmp.getWindowsInformation());
		}
		return false;
	}

	@Override
	public int compareTo(WindowsBasedBasicWritable tp) {
		return  windowsInfo.compareTo(tp.getWindows());
	}
}
