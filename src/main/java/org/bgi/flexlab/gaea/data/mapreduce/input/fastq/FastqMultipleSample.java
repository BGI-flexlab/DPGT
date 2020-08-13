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
package org.bgi.flexlab.gaea.data.mapreduce.input.fastq;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.LineReader;

import java.io.IOException;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class FastqMultipleSample {
	private Map<String,FastqSample> sampleList = new ConcurrentHashMap<String,FastqSample>();
	private boolean sequenceType = false;
	
	private int sampleNumber;
	private int mismatch = 2;
	
	public FastqMultipleSample(String list, boolean keepFile) throws IOException {
		Configuration conf = new Configuration();
		Path file = new Path(list);
		FileSystem fs = file.getFileSystem(conf);
		FSDataInputStream fsdata = fs.open(file);
		LineReader reader = new LineReader(fsdata);
		
		Text text = new Text();
		while(reader.readLine(text) != 0) {
			String line = text.toString();
			if(line.length() == 0) {
				continue;
			}
			
			FastqSample slist = new FastqSample();
			
			if(slist.setSampleList(line, keepFile)) {
				if(slist.getIndex() != null){
					sequenceType = true;
					if(slist.getFastq1() != null)
							sampleList.put(slist.getIndex(), slist);
					if(slist.getFastq2() != null)
							sampleList.put(slist.getIndex(), slist);
				}
				else
					sampleList.put(String.valueOf(slist.getId()), slist);
				sampleNumber++;
			}
		}
		reader.close();
	}
	
	public FastqSample getID(String fqName) {
		for(FastqSample slist : sampleList.values()) {
			if((slist.getFastq1() != null && fqName.contains(slist.getFastq1())) || (slist.getFastq2() != null && fqName.contains(slist.getFastq2()))) {
				return slist;
			}
		}
		return null;
	}
	
	public FastqSample getID(String index,boolean check){
		if(sampleList.containsKey(index))
			return sampleList.get(index);
		else{
			String indexFind = null;
			int minMismatch = Integer.MAX_VALUE,len = index.length(),tmpMiss,i;
			for(String key : sampleList.keySet()){
				tmpMiss = 0;
				for(i=0;i<len;i++){
					if(index.charAt(i) != key.charAt(i))
						tmpMiss++;
					if(tmpMiss > mismatch)
						break;
				}
				if(tmpMiss <= mismatch){
					if(tmpMiss < minMismatch){
						indexFind = key;
						minMismatch = tmpMiss;
					}
				}
			}
			if(minMismatch == Integer.MAX_VALUE)
				return null;
			return sampleList.get(indexFind);
		}
	}
	
	public void setMismatch(int miss){
		mismatch = miss;
	}
	
	public String getFileNameForId(int id) {
		String index = String.valueOf(id);
		if(sampleList.containsKey(index)) 
			return sampleList.get(index).getFileName();
		else
			return null;
	}

	/**
	 * @return the smapleNum
	 */
	public int getSampleNumber() {
		return sampleNumber;
	}
	
	public Map<String, FastqSample> getSampleList() {
		return sampleList;
	}
	
	public boolean getSequenceType(){
		return sequenceType;
	}
}
