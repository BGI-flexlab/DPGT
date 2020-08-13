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
package org.bgi.flexlab.gaea.util;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.mapreduce.util.HdfsFileManager;

import java.io.IOException;

public class FileIterator {
	private Configuration conf=new Configuration();
	private FileSystem fs;
	private LineReader reader=null;
	private Text value=null;
	private Path path;
	public FileIterator(String file) throws IOException {
		this.path=new Path(file);
		read();
	}
	public FileIterator(Path p) throws IOException {
		this.path=p;
		read();
	}
	
	@SuppressWarnings("deprecation")
	private void read() throws IOException {
		fs=HdfsFileManager.getFileSystem(path, conf);
		FSDataInputStream in;
		if(fs.getLength(path)==0)
			return;
		in = fs.open(path);
		reader = new LineReader(in, conf);	
	}
	
	public boolean hasNext() 
	{
		if(value!=null)
			return true;
		if(reader==null)
			return false;
		value=new Text();
		try {
			if(reader.readLine(value) > 0 && value.getLength() != 0) {
				return true;
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return false;
	}
	public Text next()
	{
		if(value==null) {
			hasNext();
		}
		Text current=value;
		value=null;
		return current;
	}
	
	
	public LineReader getReader() {
		return reader;
	}
	
	public void close() {
		if(reader!=null) {
			try {
				reader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
}
