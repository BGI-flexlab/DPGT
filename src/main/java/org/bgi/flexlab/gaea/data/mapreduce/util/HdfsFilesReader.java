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
package org.bgi.flexlab.gaea.data.mapreduce.util;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.exception.FileNotExistException;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader.HeaderPathFilter;
import org.bgi.flexlab.gaea.util.GaeaFilesReader;

import java.io.IOException;
import java.util.ArrayList;

public class HdfsFilesReader extends GaeaFilesReader{
	private Configuration conf;
	private LineReader lineReader = null;
	private ArrayList<Path> files = null;
	private FileSystem fs = null;
	
	public HdfsFilesReader(){
		this(new Configuration());
	}
	
	public HdfsFilesReader(Configuration conf){
		this.conf = conf;
		files = new ArrayList<Path>();
		currentFileIndex = 0;
		currentLine = null;
	}

	@Override
	public void traversal(String path){
		traversal(path,new HeaderPathFilter());
	}
	
	public void traversal(String path,PathFilter pathFilter) {
		Path p = new Path(path);
		
		fs = HdfsFileManager.getFileSystem(p, conf);
		FileStatus status = null;
		try {
			status = fs.getFileStatus(p);
		} catch (IOException e2) {
			throw new FileNotExistException(p.getName());
		}
		
		if(status.isFile()){
			if(!filter(p.getName()))
				files.add(p);
		}else{
			FileStatus[] stats = null;
			try{
				stats = fs.listStatus(p,pathFilter);
			}catch(IOException e){
				throw new RuntimeException(e.toString());
			}

			for (FileStatus file : stats) {
				if(!file.isFile()){
					traversal(file.toString(),pathFilter);
				}else{
					if(!filter(file.getPath().getName()))
						files.add(file.getPath());
				}
			}
		}
		
		if(size() == 0)
			return;
		FSDataInputStream currInput;
		try {
			currInput = fs.open(files.get(0));
			lineReader = new LineReader(currInput,conf);
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}

	@Override
	public boolean hasNext() {
		if(lineReader != null){
			Text line = new Text();
			try {
				if(lineReader.readLine(line) > 0){
					currentLine = line.toString();
					return true;
				}else{
					lineReader.close();
					currentFileIndex++;
					
					if(currentFileIndex < size()){
						FSDataInputStream currInput = fs.open(files.get(currentFileIndex));
						lineReader = new LineReader(currInput,conf);
						if(lineReader.readLine(line) > 0){
							currentLine = line.toString();
							return true;
						}
					}
				}
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}
		}
		
		currentLine = null;
		return false;
	}
	
	public void delete(){
		for(Path path : files){
			fs = HdfsFileManager.getFileSystem(path, conf);
			try {
				if(fs.exists(path)){
					fs.delete(path);
				}
			} catch (IOException e1) {
				throw new RuntimeException(e1.toString());
			}
		}
	}

	@Override
	public void clear() {
		if(lineReader != null){
			try {
				lineReader.close();
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}
		}
		files.clear();
		currentLine = null;
		currentFileIndex = 0;
	}

	@Override
	protected int size() {
		return files.size();
	}
}
