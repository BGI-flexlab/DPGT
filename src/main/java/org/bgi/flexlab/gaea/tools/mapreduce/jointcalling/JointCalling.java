package org.bgi.flexlab.gaea.tools.mapreduce.jointcalling;

import java.io.*;
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
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.GzipCodec;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.NLineInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.VCFHdfsWriter;
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

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableFileStream;

import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;

public class JointCalling extends ToolsRunner{
	
	public final static String INPUT_ORDER = "input.name.order";
	public final static String INPUT_LIST="input.gvcf.list";//added by gc
	public final static String Window_File="window.file.path";//window file path 
	public final static String Real_Mappers="mapper.number";
	private String HEADER_DEFAULT_PATH = "vcfheader";
	
	private String MERGER_HEADER_INFO = "vcfheaderinfo";
	public OutputStream win_out=null;
	private VCFHeader header = null;
	private LinkedHashMap<Integer, String> contigs = null;
	private boolean mapLargeMode=false;
	private boolean mapLargeMode2=false;
	private boolean mapLargeMode3=false;
	private boolean mapLargeMode4=false;
	public JointCalling(){
		this.toolsDescription = "joing calling for gvcfs";
	}
	
	private Set<String> getSampleList(Set<VCFHeader> headers){
		Set<String> samples = new TreeSet<String>();
		for(VCFHeader header : headers){
			for ( String sample : header.getGenotypeSamples() ) {
				samples.add(GaeaGvcfVariantContextUtils.mergedSampleName(null, sample, false));
			}
		}
		
		return samples;
	}
	
	private VCFHeader getVCFHeaderFromInput(Set<VCFHeader> headers) throws IOException{
        Set<String> samples = getSampleList(headers);
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(headers, true);
        VCFHeader vcfHeader = new VCFHeader(headerLines, samples);
        
        headers.clear();
        samples.clear();
        headerLines.clear();
        
        return vcfHeader;
	}

	@Override
	public int run(String[] args) throws Exception {
		BioJob job = BioJob.getInstance();
        Configuration conf = job.getConfiguration();
//        conf.setBoolean(Job.MAP_OUTPUT_COMPRESS, true);
//        conf.setClass(Job.MAP_OUTPUT_COMPRESS_CODEC, GzipCodec.class, CompressionCodec.class);
        String[] remainArgs = remainArgs(args, conf);
        JointCallingOptions options = new JointCallingOptions();
        options.parse(remainArgs);//获得header等信息
        options.setHadoopConf(remainArgs, conf);
        mapLargeMode=options.getMapperMode();
        conf.set(KeyIgnoringVCFOutputFormat.OUTPUT_VCF_FORMAT_PROPERTY, options.getOuptputFormat().toString());
        conf.setBoolean(GaeaVCFOutputFormat.HEADER_MODIFY, true);
        conf.set("mapreduce.reduce.shuffle.maxfetchfailures", "30");
        Logger logger = LoggerFactory.getLogger(VCFHeaderReader.class);
        MultipleVCFHeaderForJointCalling multiVcfHeader = new MultipleVCFHeaderForJointCalling();
        logger.warn("before get Header");
        if(options.getVcfHeaderFile() != null) {
//            conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVcfHeaderFile());
//            multiVcfHeader.headersConfig(new Path(options.getVcfHeaderFile()), options.getVCFHeaderOutput()+"/vcfHeaders", conf);
        }else {
			FileSystem fs=FileSystem.get(conf);
			Path vcfheaderPath=new Path(options.getVCFHeaderOutput()+"/"+HEADER_DEFAULT_PATH);
			Path vcfheaderinfoPath=new Path(options.getVCFHeaderOutput()+"/"+MERGER_HEADER_INFO);
			if(fs.exists(vcfheaderPath) && fs.exists(vcfheaderinfoPath)){
				System.out.println("vcfheader exists");
				conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVCFHeaderOutput() + "/" + HEADER_DEFAULT_PATH);
				conf.set(MERGER_HEADER_INFO, options.getVCFHeaderOutput() + "/" + MERGER_HEADER_INFO);
			}else{
				multiVcfHeader.headersConfig(options.getInput(), options.getVCFHeaderOutput(), conf);
			}
//        	conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVCFHeaderOutput() + "/"+HEADER_DEFAULT_PATH);
//			conf.set(MERGER_HEADER_INFO, options.getVCFHeaderOutpusampleNamePatht()+"/"+MERGER_HEADER_INFO);
            //multiVcfHeader.headersConfig(options.getInput(), options.getVCFHeaderOutput(), conf);
        }

        Path inputList=new Path(options.getInputList());
        ArrayList<String> inputListArr=new ArrayList<>();
        FileSystem inputListFs=inputList.getFileSystem(conf);
        FSDataInputStream  inputListReader=inputListFs.open(inputList);
        String tmpLine;
        Map<String,String> pathSample=new HashMap<>();
        Path vcfHeaderInfo=new Path(conf.get(MERGER_HEADER_INFO));
        FileSystem vhiFs=vcfHeaderInfo.getFileSystem(conf);
        BufferedReader indexReader = new BufferedReader(new InputStreamReader(vhiFs.open(vcfHeaderInfo)));
        while((tmpLine=indexReader.readLine())!=null) {
        	String[] eles=tmpLine.split("\t");
			if(eles.length!=3) {
				logger.error("vcfheaderinfo file format error");
			}
			String name;
			if(eles[2].endsWith(",")) {
				name=eles[2].substring(0,eles[2].length()-1);
			}else {
				name=eles[2];
			}
			pathSample.put(eles[0], name);
        }
        indexReader.close();
        HashMap<String,String> relative2Absolute=new HashMap<>();
        while((tmpLine=inputListReader.readLine())!=null) {
        	String[] eles=tmpLine.split("/");
        	String rName=eles[eles.length-1];
        	String sampleName=pathSample.get(rName);

			relative2Absolute.put(rName,tmpLine);
        }
        inputListReader.close();
        logger.warn("after get Header");

        //conf.set("mapreduce.reduce.shuffle.memory.limit.percent", "0.1");
        conf.set("mapreduce.reduce.shuffle.input.buffer.percent","0.1");
        conf.set(INPUT_LIST,options.getInputList()); 
        //create windows bed file based on window size
        conf.set(Window_File, options.getTmpOut()+"/windows.bed");

        File tmpOutDir=new File(options.getTmpOut());
        if(!tmpOutDir.exists()) {
        	tmpOutDir.mkdirs();
        }
        String win_out_file=options.getTmpOut()+"/windows.bed";
        //Path raw_win_file=new Path(options.getWinFile());//get windows file path from command
        if(win_out_file.startsWith("file://")) {
        	win_out_file=win_out_file.substring(7);
        }
        File winOutFile=new File(win_out_file);
        if(!winOutFile.exists()) {
        	winOutFile.createNewFile();
        }
        FileWriter win_out=new FileWriter(winOutFile); //output stream ready for write
        Path path = new Path(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP));
		SeekableStream in = WrapSeekable.openPath(path.getFileSystem(conf), path);
		//SeekableStream in=new SeekableFileStream(new File(options.getTmpOut()+"/virtual.vcf"));
		logger.warn("before readHeader");
		header = VCFHeaderReader.readHeaderFrom(in);
		logger.warn("after readHeader");
		in.close();
		if(header == null)
			throw new RuntimeException("header is null !!!");
        contigs = new LinkedHashMap<>();
        for (VCFContigHeaderLine line : header.getContigLines()) {
			contigs.put(line.getContigIndex(), line.getID());
		}
        int window_size=options.getWindowsSize();
        int total_bytes=0;
        for(Map.Entry<Integer, String> entry:contigs.entrySet()) {
        	//System.out.println(entry.getKey()+"\t"+entry.getValue());
        	String chr=entry.getValue();
        	int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
        	int start=1;
        	int end=-1;
        	while(true) {
        		end=start+window_size-1;
        		if(end>contigLength) {
        			end=contigLength;
        		}
        		String write_line=chr+"\t"+start+"\t"+end+"\n";
        		total_bytes+=write_line.length();
        		win_out.write(write_line);
        		start+=window_size;
        		if(start>contigLength) {
        			break;
        		}
        	}
        }
        win_out.close();
		Path vcfheaderinfoPath = new Path(conf.get(MERGER_HEADER_INFO));
		BufferedReader headerInfoReader=new BufferedReader(new InputStreamReader(vcfheaderinfoPath.getFileSystem(conf).open(vcfheaderinfoPath)));
		String headerInfoLine;
		HashMap<String,String> sampleNamePath=new HashMap<>();
		while((headerInfoLine=headerInfoReader.readLine())!=null){
			String[] eles=headerInfoLine.split("\t");
			if(eles[2].endsWith(",")){
				eles[2]=eles[2].substring(0,eles[2].length()-1);
			}
			sampleNamePath.put(eles[2],eles[0]);
//                vcfHeaderWriter.write(eles[0]+"\t"+inputIndex+"\t"+eles[2]+"\n");
//                sampleIndex.put(eles[2],inputIndex);
//                pathSample.put(eles[0], eles[2]);
//                inputIndex++;
		}
		headerInfoReader.close();

		vcfheaderinfoPath.getFileSystem(conf).delete(vcfheaderinfoPath);
		File rawInput=new File(options.getInputList());
		String sortedInputGvcfList=rawInput.getParent()+"/sorted."+rawInput.getName();
		if(sortedInputGvcfList.startsWith("file:")){
			sortedInputGvcfList=sortedInputGvcfList.substring(5);
		}

		BufferedWriter vcfInfoWriter=new BufferedWriter(new OutputStreamWriter(vcfheaderinfoPath.getFileSystem(conf).create(vcfheaderinfoPath)));
		BufferedWriter sortedInputWriter=new BufferedWriter(new FileWriter(sortedInputGvcfList));
		int inputIndex=0;
		for(String sampleName:header.getSampleNamesInOrder()){
			if(!sampleNamePath.containsKey(sampleName)){
				logger.error("code error");
				System.exit(1);
			}
			inputListArr.add(sampleName);
			sortedInputWriter.write(relative2Absolute.get(sampleNamePath.get(sampleName))+"\n");
			vcfInfoWriter.write(sampleNamePath.get(sampleName)+"\t"+inputIndex+"\t"+sampleName+"\n");
			inputIndex++;
		}
		vcfInfoWriter.close();
		sortedInputWriter.close();
		sortedInputGvcfList="file://"+sortedInputGvcfList;
		String[] inputListStringArr=new String[inputListArr.size()];
		for(int i=0;i<inputListArr.size();i++) {
			inputListStringArr[i]=inputListArr.get(i);
		}
		conf.set(INPUT_ORDER, Utils.join(",", inputListStringArr));
        int split_bytes=total_bytes/options.getMapperNumber();
        if(split_bytes<50) {
        	split_bytes=50;
        }
        int samples_size=0;
        total_bytes=0;
        BufferedReader gvcflist_reader=new BufferedReader(new FileReader(options.getInputList().substring(7)));
        String gvcf_path=null;
        while((gvcf_path=gvcflist_reader.readLine())!=null) {
        	samples_size++;
        	total_bytes+=gvcf_path.length()+7;
        }
        gvcflist_reader.close();
        split_bytes=(int)(total_bytes/options.getMapperNumber())-1;
        //conf.setLong("mapreduce.input.fileinputformat.split.maxsize", split_bytes);  //split window file by some bytes, this is a magic number which would determin the number of map task
        //System.exit(-1);
        
        //String chr = contigs.get(key.getChromosomeIndex());
		//int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
        job.setJobName("Gaea joint calling");
        //TODO:check input chr whether consistance with REF
        job.setJarByClass(JointCalling.class);
        int mapperLine=samples_size/options.getMapperNumber()<2?2:samples_size/options.getMapperNumber();
        if(samples_size%mapperLine==0) {
        	int realMappers=samples_size/mapperLine;
        	conf.set(Real_Mappers,String.valueOf(realMappers));
        }else {
        	int realMappers=(int)samples_size/mapperLine+1;
        	conf.set(Real_Mappers,String.valueOf(realMappers));
        }
		if(mapLargeMode) {
			job.setWindowsBasicMapperClass(JointCallingMapperHuge.class, options.getWindowsSize(),0);
//        				job.setMapperClass(JointCallingMapperHuge.class);
			FileInputFormat.setInputPaths(job,new Path(sortedInputGvcfList));
			job.setInputFormatClass(NLineInputFormat.class);
			NLineInputFormat.setNumLinesPerSplit(job, mapperLine);
		}else {
			FileInputFormat.setInputPaths(job, options.getInput().toArray(new Path[options.getInput().size()]));
			job.setWindowsBasicMapperClass(JointCallingMapper.class, options.getWindowsSize(),0);
			job.setInputFormatClass(JointCallingVCFInputFormat.class);
		}
        job.setReducerClass(JointCallingReducer.class);
        job.setNumReduceTasks(options.getReducerNumber());
        job.setOutputKeyValue(WindowsBasedWritable.class,VariantContextWritable.class, NullWritable.class, VariantContextWritable.class);
		job.setOutputFormatClass(GaeaVCFOutputFormat.class);
		FileOutputFormat.setOutputPath(job, new Path(options.getOutput()));
		
		//job.setOutputFormatClass(GaeaVCFOutputFormat.class);
		return job.waitForCompletion(true) ? 0 : 1;
	}
	public ArrayList<Integer> getMeanArr(int samples,int mappers){
		ArrayList<Integer> arr=new ArrayList<Integer>();
		int mean1=0;
		int mean2=0;
		if(samples<=mappers*2) {
			mean1=2;
		}else {
			mean1=samples/mappers;
		}
		if(samples % mappers==0) {
			mean2=mean1;
		}else {
			mean2=mean1+1;
		}
		int mean2_num=samples-mappers*mean1;
		int mean1_num=mean2*mappers-samples;
		for(int i=0;i<mean1_num;i++) {
			arr.add(mean1);
		}
		for(int i=0;i<mean2_num;i++) {
			arr.add(mean2);
		}
		return arr;
	}
}
