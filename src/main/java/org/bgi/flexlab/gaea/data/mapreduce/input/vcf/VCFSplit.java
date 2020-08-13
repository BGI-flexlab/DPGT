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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2010 Aalto University 
 *
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.mapreduce.input.vcf;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

public class VCFSplit {

	private SAMFileHeader mFileHeader;
	private String output;
	private FileSystem fs;
	private Configuration conf=new Configuration();
	private HashMap<String,FSDataOutputStream> sample;
	private FSDataOutputStream currentOutput;
	private boolean headerHasWrite=false;
	
	
	private final String VCFHeaderEndLineTag="#CHROM";
	private final String VCFHeaderStartTag="#";
	private final String SampleTag="SAMPLENAME:";
	public VCFSplit(String output, Configuration conf) {
		this.output=output;
		try {
			mFileHeader=SamHdfsFileHeader.getHeader(conf);
		} catch (Exception e) {
			e.printStackTrace();
		}
		sample=new HashMap<String,FSDataOutputStream>();
		Set<String> sampleName=new HashSet<String>();
		List<SAMReadGroupRecord> rg= mFileHeader.getReadGroups();
		for(SAMReadGroupRecord r:rg) {
			sampleName.add(r.getSample());
		}
		fs=getFileSystem(output,conf);
		for(String s:sampleName) {
			FSDataOutputStream data = null;
			try {
				data = fs.create(new Path(getOutputName(s)));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			sample.put(s,data);
		}
	}
	
	public void split() throws IOException {
		Path p=new Path(output+"/vcf");
		read(p);
	}
	
	private void read(Path p) throws IOException {
		if(p.getName().startsWith("_")) {
			return;
		}
		if(!fs.exists(p)) {
			throw new IOException("the input path is no exist!");
		}
		if(fs.isFile(p)) {
			readFile(p);
		} else {
			FileStatus stats[]=fs.listStatus(p);
			for(FileStatus state:stats) { 
				Path statePath=state.getPath();
				if(statePath.getName().startsWith("_")) {
					continue;
				}
				if(fs.isFile(statePath)) {
					readFile(statePath);
				} else {
					read(statePath);
				}
			}
		}
	}
	
	private void readFile(Path p) {
		StringBuilder sb=null;
		if(!headerHasWrite) {
			sb=new StringBuilder();
		}
		try {
			FSDataInputStream table=fs.open(p);
			LineReader lineReader = new LineReader(table, conf);
			Text line = new Text();
			String tempString=null;
			while(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				tempString=line.toString();
				if(tempString.startsWith(VCFHeaderStartTag)) {
					if(headerHasWrite) {
						continue;
					}
					sb.append(tempString.trim());
					sb.append("\n");
					if(tempString.startsWith(VCFHeaderEndLineTag)) {
						writeHeader(sb.toString().trim());
						headerHasWrite=true;
					}
				} else if (tempString.startsWith(SampleTag)) {
					String sampleName=tempString.split(":")[1];
					currentOutput=sample.get(sampleName);
				} else {
					write("\n");
					if(tempString.startsWith("chr1\t179462149")) {
						System.out.println("debug:"+tempString);
					}
					write(tempString);
				}
			}
			lineReader.close();
			table.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private void write(String value) {
		try {
			currentOutput.write(value.getBytes());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void writeHeader(String header) {
		byte[] h=header.getBytes();
		for(String s:sample.keySet()) {
			try {
				header = header + "\t" + s;
				h=header.getBytes();
				sample.get(s).write(h);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	private String getOutputName(String sample) {
		if(this.output.endsWith("/")) {
			return this.output+sample+".vcf";
		} else {
			return this.output+"/"+sample+".vcf";
		}
	}
	
	private static FileSystem getFileSystem(String path,Configuration conf) {
 		FileSystem fs=null;
 		try {
 			Path fsPath = new Path(path);
 			fs = fsPath.getFileSystem(conf);
		} catch (Exception e) {
			// TODO: handle exception
		}
		return fs;
 	}
	
	public void close() {
		if(sample==null) {
			return;
		}
		for(Entry<String, FSDataOutputStream> out:sample.entrySet()) {
			try {
				out.getValue().close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		sample=null;
	}
}
