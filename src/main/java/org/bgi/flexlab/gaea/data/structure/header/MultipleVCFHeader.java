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
package org.bgi.flexlab.gaea.data.structure.header;

import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.Job;

import java.io.IOException;
import java.io.Serializable;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class MultipleVCFHeader extends  GaeaVCFHeader implements Serializable{
	/**
	 * serial ID
	 */
	private static final long serialVersionUID = -5677604795673775528L;

	/**
	 * fileName 2 ID
	 */
	private Map<String, Integer> fileName2ID = new ConcurrentHashMap<String, Integer>();
	
	/**
	 * ID to header
	 */
	private Map<Integer, SingleVCFHeader> ID2SingleVcfHeader = new ConcurrentHashMap<Integer, SingleVCFHeader>();
	
	/**
	 * global ID
	 */
	private int id = 0;
	
	/**
	 * get vcf header
	 * @param id
	 * @return header String
	 */
	public VCFHeader getVcfHeader(int id) {
		return ID2SingleVcfHeader.get(id).getHeader();
	}
	
	/**
	 * get vcf header lines
	 * @param id
	 * @return
	 */
	public ArrayList<String> getVcfHeaderLines(int id) {
		ArrayList<String> headerLines = new ArrayList<String>();
		for(String line : ID2SingleVcfHeader.get(id).getHeaderInfoStringLines(null)){
			headerLines.add(line);
		}
		return headerLines;
	}
	
	/**
	 * get sample number of id file.
	 * @param id
	 * @return
	 */
	public int getSampleNum(int id) {
		return ID2SingleVcfHeader.get(id).getSampleNames().size();
	}
	
	/**
	 * get samples string of id file
	 * @param id
	 * @return
	 */
	public List<String> getSampleNames(int id) {
		return ID2SingleVcfHeader.get(id).getSampleNames();
	}
	
	/**
	 * get id from fileName, this function is for multi-vcf reader
	 * @param filePathName
	 * @return
	 */
	public int getId(String filePathName) {
		//filePathName = formatFilePath(filePathName);
		//System.err.println("getID:" + filePathName);
		if(fileName2ID.containsKey(filePathName)) {
			return fileName2ID.get(filePathName);
		} else {
			throw new RuntimeException("this file is not in inputs!");
		}
	}
	
	public String getFile(int id) {
		for(String file : fileName2ID.keySet()) {
			if(id == fileName2ID.get(file)){
				return file;
			}
		}
		throw new RuntimeException("no such id in VCFHeader!");
	}
	
	/**
	 * read single vcf file
	 * @param vcf
	 * @param conf
	 * @throws IOException
	 */
	private void readVcfHeader(Path vcf, Configuration conf) throws IOException {
		SingleVCFHeader singleVcfHeader = new SingleVCFHeader();
		singleVcfHeader.parseHeader(vcf, null, conf);
		ID2SingleVcfHeader.put(id, singleVcfHeader);
		
		//String filePathName = formatFilePath(vcf.toString());
		fileName2ID.put(vcf.toString(), id);
		
		id++;
	}
	
	@SuppressWarnings("unused")
	private String formatFilePath(String filePathName) {
		if(filePathName.startsWith("file:///")) {
			filePathName = filePathName.substring(7);
		} else {
			if(filePathName.startsWith("file:/")) {
				filePathName = filePathName.substring(5);
			}
		}
		return filePathName.trim();
	}
	
	public void mergeHeader(Path inputPath, String output, Job job, boolean distributeCacheHeader) {
		Configuration conf = job.getConfiguration();
		try {
			FileSystem fs = inputPath.getFileSystem(conf);
			fs = inputPath.getFileSystem(conf);
			if (!fs.exists(inputPath)) {
				System.out.println("Input File Path is not exist! Please check input var.");
				System.exit(-1);
			}
			if (fs.isFile(inputPath)) {
				if(validPath(inputPath, fs)){
					readVcfHeader(inputPath, conf);	
				}	
			}else {		
				FileStatus stats[]=fs.listStatus(inputPath);
				
				for (FileStatus file : stats) {
					Path filePath = file.getPath();
					mergeHeader(filePath, output, job, distributeCacheHeader);
				}
			}
			fs.close();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		if(distributeCacheHeader){
			distributeCacheVcfHeader(output, job, conf);
		} else {
			writeHeaderToHDFS(output,conf);
		}
	}
	
	private boolean validPath(Path inputPath, FileSystem fs) throws IOException{
		return (!inputPath.getName().startsWith("_")) && (fs.getFileStatus(inputPath).getLen() != 0);
	}
	
	public boolean distributeCacheVcfHeader(String outputPath, Job job, Configuration conf) {
		writeHeaderToHDFS(outputPath, conf);
		try {
			job.addCacheFile(new URI(conf.get(GaeaVCFHeader.VCF_HEADER_PROPERTY) + "#VcfHeaderObj"));
		} catch (URISyntaxException e) {
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	public int getFileNum() {
		return fileName2ID.size();
	}

	public Map<String, Integer> getFileName2ID() {
		return fileName2ID;
	}

	public Map<Integer, SingleVCFHeader> getID2SingleVcfHeader() {
		return ID2SingleVcfHeader;
	}
	
}
