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

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeader.HEADER_FIELDS;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class SingleVCFHeader extends GaeaVCFHeader implements Serializable{
	/**
	 * serial id
	 */
	private static final long serialVersionUID = -5010893132552487497L;
	
	/**
	 * sample names
	 */
	private List<String> sampleNames = new ArrayList<String>();
	
	/**
	 * basic header string without sample names 
	 */
	private String headerInfo;
	
	private VCFHeader vcfHeader;
	
	/**
	 * parse header from given VCF files
	 * @param inVcf
	 * @param output
	 * @param conf
	 * @throws IOException
	 */
	public void parseHeader(Path vcf, String output, Configuration conf) throws IOException {
		readSingleHeader(vcf, conf);
		
		if(output != null) {//multi-sample do not write to HDFS, top class will handle this.
			writeHeaderToHDFS(output, conf);
		}
	}
	
	public void readHeaderFrom(Path path, FileSystem fs) throws IOException {
		SeekableStream i = WrapSeekable.openPath(fs, path);
		readHeaderFrom(i);
		i.close();
	}
	
	public void readHeaderFrom(SeekableStream in) throws IOException {
		this.setHeader(VCFHeaderReader.readHeaderFrom(in));
	}
	
	public void readSingleHeader(Path vcfPath, Configuration conf) throws IOException {
		FileSystem fs = vcfPath.getFileSystem(conf);
		if(!fs.exists(vcfPath))
			throw new RuntimeException(vcfPath.toString() + " don't exists.");
		if(!fs.isFile(vcfPath)) {
			throw new RuntimeException(vcfPath.toString() + " is not a file. GaeaSingleVcfHeader parser only support one vcf file.");
		}
		FSDataInputStream in = fs.open(vcfPath);
		AsciiLineReaderIterator it = new AsciiLineReaderIterator(new AsciiLineReader(in));
	    VCFCodec codec = new VCFCodec();
	    Object header = codec.readHeader(it);
		vcfHeader = (VCFHeader)(((FeatureCodecHeader)header).getHeaderValue());
		sampleNames.addAll(vcfHeader.getGenotypeSamples());
		buildHeaderInfo();
		it.close();
	}
	
	public void buildHeaderInfo() {
		StringBuilder headerInfo = new StringBuilder();
		for(VCFHeaderLine line : vcfHeader.getMetaDataInInputOrder()) {
			headerInfo.append(line.toString());
			headerInfo.append("\n");
		}
		headerInfo.append(VCFHeader.HEADER_INDICATOR);
		for(HEADER_FIELDS field : vcfHeader.getHeaderFields()) {
			headerInfo.append(field.toString());
			headerInfo.append("\t");
		}
		this.headerInfo = headerInfo.toString();
	}
	
	private void collectMetaInfo(StringBuilder header, String line){
		header.append(line.trim());
		header.append("\n");
	}
	
	private void collectHeader(StringBuilder header, String line){
		String[] tags = line.trim().split("\t");
		header.append(tags[0]);
		for(int i = 1; i < 9; i++) {
			header.append("\t");
			header.append(tags[i]);
		}
	}
	
	private void collectSamples(String line){
		String[] lineSplits = line.split("\t");
		for(int i = 9; i < lineSplits.length; i++) {
			sampleNames.add(lineSplits[i]);
		}
	}
	
	public String getHeaderInfoWithSample(List<String> samples) {
		if(samples == null) {
			samples = this.sampleNames;
		}
		StringBuilder sb = new StringBuilder();
		sb.append(headerInfo);
		for(String sample : samples) {
			sb.append("\t");
			sb.append(sample);
		}
		sb.append("\n");
		return sb.toString();
	}
	
	public String[] getHeaderInfoStringLines(List<String> samples) {
		return getHeaderInfoWithSample(samples).split("\n");
	}
	
	public String getSampleNames(int index) {
		return sampleNames.get(index);
	} 
	
	public List<String> getSampleNames() {
		return sampleNames;
	}
	
	public String getHeaderInfo() {
		return headerInfo;
	}

	@Deprecated
	public VCFHeaderVersion getVCFVersion(VCFHeader vcfHeader) {
		String versionLine = null;
		Set<VCFHeaderLine> vcfHeaderLineSet = vcfHeader.getMetaDataInInputOrder();
		for (VCFHeaderLine vcfHeaderLine : vcfHeaderLineSet) {
    		if(VCFHeaderVersion.isFormatString(vcfHeaderLine.getKey())){
    			versionLine = vcfHeaderLine.toString();
    			break;
    		}
		}
		return VCFHeaderVersion.getHeaderVersion("##"+versionLine);
	}

	public static VCFHeaderVersion getVCFHeaderVersion(VCFHeader vcfHeader) {
		String versionLine = null;
		Set<VCFHeaderLine> vcfHeaderLineSet = vcfHeader.getMetaDataInInputOrder();
		for (VCFHeaderLine vcfHeaderLine : vcfHeaderLineSet) {
			if(VCFHeaderVersion.isFormatString(vcfHeaderLine.getKey())){
				versionLine = vcfHeaderLine.toString();
				break;
			}
		}
		return VCFHeaderVersion.getHeaderVersion("##"+versionLine);
	}

	public VCFHeader getHeader() {
		return vcfHeader;
	}

	public void setHeader(VCFHeader header) {
		this.vcfHeader = header;
	}

}
