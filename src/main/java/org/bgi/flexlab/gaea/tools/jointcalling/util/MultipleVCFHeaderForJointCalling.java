package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.Path;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.input.vcf.VCFHdfsLoader;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.VCFHdfsWriter;
import org.bgi.flexlab.gaea.data.structure.header.GaeaVCFHeader;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFFileWriter;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalWriter;
import org.bgi.flexlab.gaea.util.FileIterator;
import org.seqdoop.hadoop_bam.LazyVCFGenotypesContext.HeaderDataCache;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.vcf.VCFHeader;

public class MultipleVCFHeaderForJointCalling extends GaeaVCFHeader implements Serializable {
	private static final long serialVersionUID = -3352974548797591555L;

	private String HEADER_DEFAULT_PATH = "vcfheader";
	
	private String MERGER_HEADER_INFO = "vcfheaderinfo";

	private int currentIndex = 0;

	//private HashMap<String, Integer> nameIndexs = new HashMap<String, Integer>();
	private HashMap<Integer, Set<String>> nameHeaders = new HashMap<Integer, Set<String>>();

	private HashSet<VCFHeader> headers = new HashSet<VCFHeader>();

	private FSDataOutputStream outputStream = null;

	private VCFHeader mergeHeader = null;
	
	//private HashMap<Integer,HeaderDataCache> vcfHeaderDateCaches = new HashMap<Integer,HeaderDataCache>();

	public MultipleVCFHeaderForJointCalling() {
	}

	/*public Set<String> getSampleList(String sampleName) {
		if (nameIndexs.containsKey(sampleName))
			return getSampleList(nameIndexs.get(sampleName));
		return null;
	}*/

	public Set<String> getSampleList(int index) {
		if (nameHeaders.containsKey(index)) {
			if(nameHeaders.get(index).size() == 0)
				return null;
			return nameHeaders.get(index);
		}
		return null;
	}

	/*
	public Integer getIndex(String name) {
		if (nameIndexs.containsKey(name))
			return nameIndexs.get(name);
		return null;
	}
	*/

	public void headersConfig(List<Path> paths, String outdir, Configuration conf) {
		conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, outdir + "/"+HEADER_DEFAULT_PATH);
		conf.set(MERGER_HEADER_INFO, outdir+"/"+MERGER_HEADER_INFO);
		
		getHeaders(paths, outdir, conf);
		writeMergeHeaders(conf);
	}
	public static Set<String> getSampleList(Set<VCFHeader> headers) {
		Set<String> samples = new TreeSet<String>();
		for (VCFHeader header : headers) {
			for (String sample : header.getGenotypeSamples()) {
				samples.add(GaeaGvcfVariantContextUtils.mergedSampleName(null, sample, false));
			}
		}

		return samples;
	}
	private VCFHeader smartMergeHeaders(Set<VCFHeader> headers) {
		Set<String> samplelists = getSampleList(headers);
		Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(headers, true);
		VCFHeader vcfHeader = new VCFHeader(headerLines, samplelists);
		headers.clear();
		samplelists.clear();
		headerLines.clear();

		return vcfHeader;
	}

	private void writeInfo(Path outputpath, Configuration conf, String fileName, StringBuilder samples) {
		if (outputStream == null) {
			try {
				outputStream = outputpath.getFileSystem(conf).create(outputpath);
			} catch (IOException e) {
				throw new UserException.CouldNotCreateOutputFile(outputpath.toString(), "create header info error.");
			}
		}

		try {
			outputStream.write((fileName + "\t" + currentIndex + "\t" + samples.toString()+"\n").getBytes());
		} catch (IOException e) {
			throw new UserException("write header info error.");
		}
	}

	public void getHeaders(List<Path> paths, String outdir, Configuration conf) {
		currentIndex = 0;
		Path outpath = new Path(outdir+"/"+MERGER_HEADER_INFO);
		for (Path p : paths) {
			VCFHdfsLoader loader = null;
			try {
				loader = new VCFHdfsLoader(p.toString());
			} catch (IllegalArgumentException | IOException e) {
				throw new UserException.BadInput("read vcf file error.");
			}

			if (loader.getHeader().getSampleNamesInOrder().size() == 0)
				throw new UserException.MalformedVCF("VCF header contains no samples!");

			StringBuilder samples = new StringBuilder();
			for (String sample : loader.getHeader().getSampleNamesInOrder()) {
				samples.append(sample + ",");
			}

			writeInfo(outpath, conf, p.getName(), samples);

			headers.add(loader.getHeader());

			if (headers.size() >= 1000) {
				VCFHeader vcfHeader = smartMergeHeaders(headers);
				headers.clear();
				headers.add(vcfHeader);
			}

			currentIndex++;
			loader.close();
		}

		if (headers.size() >= 1) {
			mergeHeader = smartMergeHeaders(headers);
			headers.clear();
		}
		
		close();
	}

	private void readHeaderInfo(String outputpath, Configuration conf) {
		try {
			FileIterator iterator = new FileIterator(outputpath);

			while (iterator.hasNext()) {
				String[] str = iterator.next().toString().split("\t");
				Set<String> samples = new HashSet<String>();
				String[] sample = str[2].split(",");
				for (String s : sample)
					samples.add(s);
				nameHeaders.put(Integer.parseInt(str[1]), samples);
			}

			iterator.close();
		} catch (IOException e) {
			throw new RuntimeException(" read header info error.");
		}
	}
	
	public int getIndex(Configuration conf,String fileName) {
		String outputpath = conf.get(MERGER_HEADER_INFO);
		int index = -1;
		try {
			FileIterator iterator = new FileIterator(outputpath);

			while (iterator.hasNext()) {
				String[] str = iterator.next().toString().split("\t");
				
				if(str[0].equals(fileName)) {
					index = Integer.parseInt(str[1]);
					break;
				}
			}

			iterator.close();
		} catch (IOException e) {
			throw new RuntimeException(" read header info error.");
		}
		
		return index;
	}

	public void readHeaders(Configuration conf) {
		String headerString = conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP);
		readHeader(new Path(headerString), conf);
		readHeaderIndex(conf);
	}
	
	public void readHeaderIndex(Configuration conf) {
		String headerInfo = conf.get(MERGER_HEADER_INFO);
		readHeaderInfo(headerInfo, conf);
		this.currentIndex = nameHeaders.size();
		//parseHeaderDataCache();
	}
	
	private void parseHeaderDataCache() {
		for(int i = 0 ; i < this.currentIndex ; i++) {
			Set<String> samples = getSampleList(i);
			if(samples == null) {
				throw new UserException("samples null!!");
			}
			
			VCFHeader header = new VCFHeader(mergeHeader.getMetaDataInInputOrder(), samples);
			
			HeaderDataCache datacache = new HeaderDataCache();
			datacache.setHeader(header);
			System.out.println(i+"\tmergeHeader:\n"+mergeHeader);
			System.out.println(i+"\theader\t"+header.getSampleNamesInOrder().toString()+"\n"+header);
			System.out.println(i+"\tdatacache:\n"+datacache);
			System.out.println(i+"\tname:\n"+datacache.toString());
			//vcfHeaderDateCaches.put(i, datacache);
			
			if(i>0) {
				System.exit(-1);
			}
		}
	}
	
//	public HashMap<Integer,HeaderDataCache> getHeaderDataCache(){
//		return this.vcfHeaderDateCaches;
//	}

	public void readHeader(Path path, Configuration conf) {
		SeekableStream in = null;
		try {
			in = WrapSeekable.openPath(path.getFileSystem(conf), path);
		} catch (IOException e1) {
			throw new UserException.CouldNotReadInputFile("cann't open path " + path.toString());
		}
		try {
			this.mergeHeader = VCFHeaderReader.readHeaderFrom(in);
			in.close();
		} catch (IOException e) {
			throw new UserException("read/close vcf file error!");
		}

	}
	
	private void writeMergeHeaders(Configuration conf) {
		VCFHdfsWriter vcfHdfsWriter = null;
		try {
			vcfHdfsWriter = new VCFHdfsWriter(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP), false, false, conf);
		} catch (IOException e) {
			throw new UserException(e.toString());
		}
		vcfHdfsWriter.writeHeader(this.mergeHeader);
		vcfHdfsWriter.close();
	}
	private void writeMergeHeadersLocal(Configuration conf) {
		VCFLocalWriter vcfHdfsWriter = null;
		try {
			vcfHdfsWriter = new VCFLocalWriter(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP), false, false);
		} catch (IOException e) {
			throw new UserException(e.toString());
		}
		vcfHdfsWriter.writeHeader(this.mergeHeader);
		vcfHdfsWriter.close();
	}
	public int getHeaderSize() {
		return this.currentIndex;
	}

	public VCFHeader getMergeHeader() {
		return mergeHeader;
	}
	public void setMergeHeader(VCFHeader mergeHeader){
		this.mergeHeader=mergeHeader;
	}
	public void setNameHeaders(HashMap<Integer, Set<String>> NameHeaders){
		this.nameHeaders.putAll(NameHeaders);
	}
	public void setCurrentIndex(int index){
		this.currentIndex=index;
	}
	public void close() {
		if (outputStream != null) {
			try {
				outputStream.close();
			} catch (IOException e) {
				throw new RuntimeException("vcf output stream close errpr. " + e.toString());
			}
		}
	}

	public void clear() {
		//nameIndexs.clear();
		nameHeaders.clear();
		headers.clear();
		//vcfHeaderDateCaches.clear();
	}

}
