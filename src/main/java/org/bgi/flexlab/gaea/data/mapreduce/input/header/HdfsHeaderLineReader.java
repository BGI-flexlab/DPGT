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
package org.bgi.flexlab.gaea.data.mapreduce.input.header;

import htsjdk.samtools.util.LineReader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;

import java.io.IOException;

public class HdfsHeaderLineReader implements LineReader{
	private int lineNumber;
	private String currentLine;
	private HdfsLineReader lineReader = null;
	
	public HdfsHeaderLineReader(Path path, Configuration conf) throws IOException{
		FileSystem fs = path.getFileSystem(conf);
		FSDataInputStream FSinput = fs.open(path);
		lineReader = new HdfsLineReader(FSinput,conf);
	}

	@Override
	public int getLineNumber() {
		return lineNumber;
	}

	@Override
	public int peek() {
		Text tmp = new Text();
		try {
			if(lineReader.readLine(tmp) == 0)
				return -1;
			else
				lineNumber++;
		} catch (IOException e) {
			e.printStackTrace();
		}
		currentLine = tmp.toString();
		return currentLine.charAt(0);
	}

	@Override
	public String readLine() {
		return currentLine;
	}

	@Override
	public void close() {
		if(lineReader == null)
			return;
		try {
			lineReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
