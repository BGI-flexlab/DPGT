package org.bgi.flexlab.gaea.tools.jointcalling.util;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.vcf.VCFHeader;
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
import org.bgi.flexlab.gaea.util.FileIterator;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;

public class MultipleVCFHeaderForJointCalling extends GaeaVCFHeader implements Serializable {
	private static final long serialVersionUID = -3352974548797591555L;

	private final String MERGER_HEADER_INFO = "vcfheaderinfo";

	private int currentIndex = 0;

	private final HashMap<Integer, Set<String>> nameHeaders = new HashMap<>();

	private final HashSet<VCFHeader> headers = new HashSet<>();

	private FSDataOutputStream outputStream = null;

	private VCFHeader mergeHeader = null;
	

	public MultipleVCFHeaderForJointCalling() {
	}


	public Set<String> getSampleList(int index) {
		if (nameHeaders.containsKey(index)) {
			if(nameHeaders.get(index).size() == 0)
				return null;
			return nameHeaders.get(index);
		}
		return null;
	}


	public void headersConfig(List<Path> paths, String outdir, Configuration conf) {
		String HEADER_DEFAULT_PATH = "vcfheader";
		conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, outdir + "/"+ HEADER_DEFAULT_PATH);
		conf.set(MERGER_HEADER_INFO, outdir+"/"+MERGER_HEADER_INFO);
		
		getHeaders(paths, outdir, conf);
		writeMergeHeaders(conf);
	}
	public static Set<String> getSampleList(Set<VCFHeader> headers) {
		Set<String> samples = new TreeSet<>();
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
			VCFHdfsLoader loader;
			try {
				loader = new VCFHdfsLoader(p.toString());
			} catch (IllegalArgumentException | IOException e) {
				throw new UserException.BadInput("read vcf file error.");
			}

			if (loader.getHeader().getSampleNamesInOrder().size() == 0)
				throw new UserException.MalformedVCF("VCF header contains no samples!");

			StringBuilder samples = new StringBuilder();
			for (String sample : loader.getHeader().getSampleNamesInOrder()) {
				samples.append(sample).append(",");
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

	private void readHeaderInfo(String outputpath) {
		try {
			FileIterator iterator = new FileIterator(outputpath);

			while (iterator.hasNext()) {
				String[] str = iterator.next().toString().split("\t");
				String[] sample = str[2].split(",");
				Set<String> samples = new HashSet<>(Arrays.asList(sample));
				nameHeaders.put(Integer.parseInt(str[1]), samples);
			}

			iterator.close();
		} catch (IOException e) {
			throw new RuntimeException(" read header info error.");
		}
	}

	public void readHeaders(Configuration conf) {
		String headerString = conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP);
		readHeader(new Path(headerString), conf);
		readHeaderIndex(conf);
	}
	
	public void readHeaderIndex(Configuration conf) {
		String headerInfo = conf.get(MERGER_HEADER_INFO);
		readHeaderInfo(headerInfo);
		this.currentIndex = nameHeaders.size();
	}



	public void readHeader(Path path, Configuration conf) {
		SeekableStream in;
		try {
			in = WrapSeekable.openPath(path.getFileSystem(conf), path);
		} catch (IOException e1) {
			throw new UserException.CouldNotReadInputFile("cann't open path " + path);
		}
		try {
			this.mergeHeader = VCFHeaderReader.readHeaderFrom(in);
			in.close();
		} catch (IOException e) {
			throw new UserException("read/close vcf file error!");
		}

	}
	
	private void writeMergeHeaders(Configuration conf) {
		VCFHdfsWriter vcfHdfsWriter;
		try {
			vcfHdfsWriter = new VCFHdfsWriter(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP), false, false, conf);
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
				throw new RuntimeException("vcf output stream close errpr. " + e);
			}
		}
	}

	public void clear() {
		nameHeaders.clear();
		headers.clear();
	}

}
