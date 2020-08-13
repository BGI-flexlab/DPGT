package org.bgi.flexlab.gaea.tools.mapreduce.jointcallingprepare;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.URI;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.GzipCodec;
import org.apache.hadoop.mapred.lib.TotalOrderPartitioner;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Partitioner;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.NLineInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.mapreduce.lib.partition.InputSampler;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.VCFHdfsWriter;
import org.bgi.flexlab.gaea.data.mapreduce.partitioner.WindowsBasedPartitioner;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.bgi.flexlab.gaea.util.FileIterator;
import org.bgi.flexlab.gaea.util.Utils;
import org.seqdoop.hadoop_bam.KeyIgnoringVCFOutputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;

public class JointCallingPrepare extends ToolsRunner {

	public final static String INPUT_ORDER = "input.name.order";
	public final static String INPUT_LIST = "input.gvcf.list";// added by gc
	public final static String Window_File = "window.file.path";// window file path
	public final static String Real_Mappers = "mapper.number";
	public final static String AccumulateIndex="refChrAccLen";
	public static Long refLength=0L;

	private String HEADER_DEFAULT_PATH = "vcfheader";

	private String MERGER_HEADER_INFO = "vcfheaderinfo";
	public OutputStream win_out = null;
	private VCFHeader header = null;
	private LinkedHashMap<Integer, String> contigs = null;
	public Long reduceBand=1L;
	public JointCallingPrepare() {
		this.toolsDescription = "Prepare for joing calling";
	}

	private Set<String> getSampleList(Set<VCFHeader> headers) {
		Set<String> samples = new TreeSet<String>();
		for (VCFHeader header : headers) {
			for (String sample : header.getGenotypeSamples()) {
				samples.add(GaeaGvcfVariantContextUtils.mergedSampleName(null, sample, false));
			}
		}

		return samples;
	}

	private VCFHeader getVCFHeaderFromInput(Set<VCFHeader> headers) throws IOException {
		Set<String> samples = getSampleList(headers);
		Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(headers, true);
		VCFHeader vcfHeader = new VCFHeader(headerLines, samples);

		headers.clear();
		samples.clear();
		headerLines.clear();

		return vcfHeader;
	}
	public class WindowsBasedGenomeAllOrderPartitioner<T> extends Partitioner<WindowsBasedWritable, T> {

		@Override
		public int getPartition(WindowsBasedWritable key, T v, int numPartitioner) {
			//whole genome
			Integer chr=key.getChromosomeIndex();
			Long pos=0L;
			for(int i=0;i<pos;i++) {
				pos+=header.getSequenceDictionary().getSequence(i).getSequenceLength();
			}
			pos+=key.getPosition().get();
			return (int)(pos / reduceBand);
		}
	}
	@Override
	public int run(String[] args) throws Exception {
		BioJob job = BioJob.getInstance();
		Configuration conf = job.getConfiguration();
		conf.setBoolean(Job.MAP_OUTPUT_COMPRESS, true);
		conf.setClass(Job.MAP_OUTPUT_COMPRESS_CODEC, GzipCodec.class, CompressionCodec.class);
		String[] remainArgs = remainArgs(args, conf);
		JointCallingPrepareOptions options = new JointCallingPrepareOptions();
		options.parse(remainArgs);// 获得header等信息
		options.setHadoopConf(remainArgs, conf);
		conf.set(KeyIgnoringVCFOutputFormat.OUTPUT_VCF_FORMAT_PROPERTY, options.getOuptputFormat().toString());
		conf.setBoolean(GaeaVCFOutputFormat.HEADER_MODIFY, true);
		conf.set("mapreduce.reduce.shuffle.maxfetchfailures", "30");
		Logger logger = LoggerFactory.getLogger(VCFHeaderReader.class);
		MultipleVCFHeaderForJointCalling multiVcfHeader = new MultipleVCFHeaderForJointCalling();
		logger.warn("before get Header");
		if (options.getVcfHeaderFile() != null) {
//            conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVcfHeaderFile());
//            multiVcfHeader.headersConfig(new Path(options.getVcfHeaderFile()), options.getVCFHeaderOutput()+"/vcfHeaders", conf);
		} else {
//            conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVCFHeaderOutput() + "/vcfFileHeader.vcf");
//			conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVCFHeaderOutput() + "/" + HEADER_DEFAULT_PATH);
//			conf.set(MERGER_HEADER_INFO, options.getVCFHeaderOutput() + "/" + MERGER_HEADER_INFO);
			FileSystem fs=FileSystem.get(conf);
			Path vcfheaderPath=new Path(options.getVCFHeaderOutput()+"/"+HEADER_DEFAULT_PATH);
			Path vcfheaderinfoPath=new Path(options.getVCFHeaderOutput()+"/"+MERGER_HEADER_INFO);
			if(fs.exists(vcfheaderPath) && fs.exists(vcfheaderinfoPath)){
				logger.warn("vcfheader exists");
				conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVCFHeaderOutput() + "/" + HEADER_DEFAULT_PATH);
				conf.set(MERGER_HEADER_INFO, options.getVCFHeaderOutput() + "/" + MERGER_HEADER_INFO);
			}else{
				multiVcfHeader.headersConfig(options.getInput(), options.getVCFHeaderOutput(), conf);
			}
			// multiVcfHeader.headersConfig(options.getInput(),
			// options.getVCFHeaderOutput(), conf);
//            VCFHeader vcfHeader = getVCFHeaderFromInput(multiVcfHeader.getHeaders());
//            VCFHdfsWriter vcfHdfsWriter = new VCFHdfsWriter(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP), false, false, conf);
//            vcfHdfsWriter.writeHeader(vcfHeader);
//            vcfHdfsWriter.close();
		}
//        String headerString = options.getOutput()+"/../vcfheaderinfo";
//		FileIterator iterator = new FileIterator(headerString);
//		System.out.println(headerString);
//		System.exit(-1);
		Path inputList = new Path(options.getInputList());
		ArrayList<String> inputListArr = new ArrayList<>();
		FileSystem inputListFs = inputList.getFileSystem(conf);
		FSDataInputStream inputListReader = inputListFs.open(inputList);
		String tmpLine;
		Map<String, String> pathSample = new HashMap<>();
		Path vcfHeaderInfo = new Path(conf.get(MERGER_HEADER_INFO));
		FileSystem vhiFs = vcfHeaderInfo.getFileSystem(conf);
		BufferedReader indexReader = new BufferedReader(new InputStreamReader(vhiFs.open(vcfHeaderInfo)));
		while ((tmpLine = indexReader.readLine()) != null) {
			String[] eles = tmpLine.split("\t");
			if (eles.length != 3) {
				logger.error("vcfheaderinfo file format error");
			}
			String name;
			if (eles[2].endsWith(",")) {
				name = eles[2].substring(0, eles[2].length() - 1);
			} else {
				name = eles[2];
			}
			pathSample.put(eles[0], name);
		}
		indexReader.close();
		while ((tmpLine = inputListReader.readLine()) != null) {
			String[] eles = tmpLine.split("/");
			String rName = eles[eles.length - 1];
			String sampleName = pathSample.get(rName);
			inputListArr.add(sampleName);
		}
		inputListReader.close();
		logger.warn("after get Header");
		String[] inputListStringArr = new String[inputListArr.size()];
		for (int i = 0; i < inputListArr.size(); i++) {
			inputListStringArr[i] = inputListArr.get(i);
		}
		conf.set(INPUT_ORDER, Utils.join(",", inputListStringArr));
		inputListArr.clear();
		// conf.set("mapreduce.reduce.shuffle.memory.limit.percent", "0.1");
		conf.set("mapreduce.reduce.shuffle.input.buffer.percent", "0.1");
		conf.set(INPUT_LIST, options.getInputList());
		// create windows bed file based on window size
		conf.set(Window_File, options.getTmpOut() + "/windows.bed");

		File tmpOutDir = new File(options.getTmpOut());
		if (!tmpOutDir.exists()) {
			tmpOutDir.mkdirs();
		}
		String win_out_file = options.getTmpOut() + "/windows.bed";
		// Path raw_win_file=new Path(options.getWinFile());//get windows file path from
		// command
		if (win_out_file.startsWith("file://")) {
			win_out_file = win_out_file.substring(7);
		}
		File winOutFile = new File(win_out_file);
		if (!winOutFile.exists()) {
			winOutFile.createNewFile();
		}
		FileWriter win_out = new FileWriter(winOutFile); // output stream ready for write
		Path path = new Path(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP));
		BufferedReader totalSampleVCF=new BufferedReader(new InputStreamReader(path.getFileSystem(conf).open(path)));
		FileWriter virtualVCF=new FileWriter(new File(options.getTmpOut()+"/virtual.vcf"));
//		Path vVCF=new Path(options.getTmpOut()+"/virtual.vcf");
//		FSDataOutputStream  virtualVCF=vVCF.getFileSystem(conf).create(vVCF);
		
		String tmpline;
		while((tmpline=totalSampleVCF.readLine())!=null) {
			if(tmpline.startsWith("#CHROM")) {
				String writeLine="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVirtualSample\n";
				virtualVCF.write(writeLine);
				break;
			}else {
				virtualVCF.write(tmpline);
				virtualVCF.write("\n");
			}
		}
		virtualVCF.close();
		SeekableStream in=new SeekableFileStream(new File(options.getTmpOut()+"/virtual.vcf"));
		header = VCFHeaderReader.readHeaderFrom(in);
		in.close();
		if (header == null)
			throw new RuntimeException("header is null !!!");
		contigs = new LinkedHashMap<>();
		for (VCFContigHeaderLine line : header.getContigLines()) {
			contigs.put(line.getContigIndex(), line.getID());
		}
		int window_size = options.getWindowsSize();
		int total_bytes = 0;
		for (Map.Entry<Integer, String> entry : contigs.entrySet()) {
			// System.out.println(entry.getKey()+"\t"+entry.getValue());
			String chr = entry.getValue();
			int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();

			refLength+=contigLength;
			int start = 1;
			int end = -1;
			while (true) {
				end = start + window_size - 1;
				if (end > contigLength) {
					end = contigLength;
				}
				String write_line = chr + "\t" + start + "\t" + end + "\n";
				total_bytes += write_line.length();
				win_out.write(write_line);
				start += window_size;
				if (start > contigLength) {
					break;
				}
			}
		}
		win_out.close();
		int split_bytes = total_bytes / options.getMapperNumber();
		if (split_bytes < 50) {
			split_bytes = 50;
		}
		int samples_size = 0;
		total_bytes = 0;
		BufferedReader gvcflist_reader = new BufferedReader(new FileReader(options.getInputList().substring(7)));
		String gvcf_path = null;
		while ((gvcf_path = gvcflist_reader.readLine()) != null) {
			samples_size++;
			total_bytes += gvcf_path.length() + 7;
		}
		gvcflist_reader.close();
		split_bytes = (int) (total_bytes / options.getMapperNumber()) - 1;
		// conf.setLong("mapreduce.input.fileinputformat.split.maxsize", split_bytes);
		// //split window file by some bytes, this is a magic number which would
		// determin the number of map task
		// System.exit(-1);

		// String chr = contigs.get(key.getChromosomeIndex());
		// int contigLength =
		// header.getSequenceDictionary().getSequence(chr).getSequenceLength();
		job.setJobName("Gaea joint calling preparation");
		// TODO:check input chr whether consistance with REF
		job.setJarByClass(JointCallingPrepare.class);
		int mapperLine = samples_size / options.getMapperNumber() < 2 ? 2 : samples_size / options.getMapperNumber();
		if (samples_size % mapperLine == 0) {
			int realMappers = samples_size / mapperLine;
			conf.set(Real_Mappers, String.valueOf(realMappers));
		} else {
			int realMappers = (int) samples_size / mapperLine + 1;
			conf.set(Real_Mappers, String.valueOf(realMappers));
		}
//		String tmpDir = options.getTmpOut() + "/TMP";
//		if (tmpDir.startsWith("file://")) {
//			tmpDir = tmpDir.substring(7);
//		}
//		File tmpDirfile = new File(tmpDir);
//		if (!tmpDirfile.exists()) {
//			tmpDirfile.mkdirs();
//		}
//		String tmpDir2 = options.getTmpOut() + "/TMPDONE";
//		if (tmpDir2.startsWith("file://")) {
//			tmpDir2 = tmpDir2.substring(7);
//		}
//		File tmpDirfile2 = new File(tmpDir2);
//		if (!tmpDirfile2.exists()) {
//			tmpDirfile2.mkdirs();
//		}
		job.setWindowsBasicMapperClass(JointCallingPrepareMapper.class, options.getWindowsSize(), 0);
		FileInputFormat.setInputPaths(job, new Path(conf.get(INPUT_LIST)));
		job.setInputFormatClass(NLineInputFormat.class);
		NLineInputFormat.setNumLinesPerSplit(job, mapperLine);
		job.setReducerClass(JointCallingPrepareReducer.class);
		job.setNumReduceTasks(options.getReducerNumber());
		reduceBand=refLength%options.getReducerNumber()==0?(Long)(refLength/options.getReducerNumber()):(Long)(refLength/(options.getReducerNumber()+1));
		//System.out.println("reduce band:\t"+reduceBand);
		//全排序
		//System.out.println(job.getSortComparator().getClass());
		job.setPartitionerClass(WindowsBasedTestPartitioner.class);
		//InputSampler.Sampler<WindowsBasedWritable, Text> sampler=new InputSampler.RandomSampler<WindowsBasedWritable,Text>(0.1, 1000, 10);
		
	//InputSampler.writePartitionFile(job, sampler);
//		String partitionFile=TotalOrderPartitioner.getPartitionFile(conf);
//		URI partitionUri=new URI(partitionFile);
//		job.addCacheFile(partitionUri);
		
		job.setOutputKeyValue(WindowsBasedWritable.class, Text.class, NullWritable.class,
				Text.class);
		job.setOutputFormatClass(TextOutputFormat.class);
		FileOutputFormat.setOutputPath(job, new Path(options.getOutput()));
		boolean returnValue=job.waitForCompletion(true);
		return  returnValue?1:0;
	}

}
