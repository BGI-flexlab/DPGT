package org.bgi.flexlab.gaea.tools.mapreduce.jointcallingprepare;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.URI;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

import javax.ws.rs.core.Variant;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.data.structure.dbsnp.DbsnpShare;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.data.structure.variant.VariantContextMerger;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.variant.filter.VariantRegionFilter;
import org.bgi.flexlab.gaea.tools.jointcalling.JointCallingEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.VariantAnnotatorEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.bgi.flexlab.gaea.tools.jointcalling.util.ReferenceConfidenceVariantContextMerger;
import org.bgi.flexlab.gaea.util.Utils;
import org.seqdoop.hadoop_bam.LazyParsingGenotypesContext;
import org.seqdoop.hadoop_bam.LazyVCFGenotypesContext;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.seqdoop.hadoop_bam.LazyVCFGenotypesContext.HeaderDataCache;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFUtils;

public class JointCallingPrepareMapper extends
	Mapper<LongWritable, Text, WindowsBasedWritable, Text>{
	public static HashMap<Integer,Long> accumulate=new HashMap<>();
	//private int windowSize = 10000;
	//private int thread_num=10;
	public static Long refLength=0L;
	private Logger logger = LoggerFactory.getLogger(JointCallingPrepareMapper.class);
	private int run_times=0;
	private int processed_samples=0;
	private HashMap<String,Integer> chrIndexs = null;
	public static LinkedHashSet<String> gvcflist=new LinkedHashSet<String>();
	private WindowsBasedWritable outKey = new WindowsBasedWritable();
	private ArrayList<VCFCodec> codec = new ArrayList<VCFCodec>();
	public static VCFHeader merged_header = null;
	private VCFHeader mapMergedHeader=null;
	private String mapSMtag=null;
	private Integer mapSMtagInt=0;
	private String mapSMfirstTag=null;
	private ArrayList<VCFHeader> sample_headers=new ArrayList<VCFHeader>();
	private ArrayList<String> sampleNames=new ArrayList<String>();
	private String gvcfList=null;
	private String tmpDir=null;
	private String tmpDir2=null;
	private String tmpOutStr=null;
	private int realMapperNumber=0;
	Configuration conf;
	private int windowSize;
	private HashMap<Integer, String> contigs = new HashMap<Integer, String>();
	private Map<String, Integer> contigDict = new HashMap<String, Integer>();
	//private JointCallingEngine engine = null;
	private ArrayList<VariantContextWritable> outValue = new ArrayList<VariantContextWritable>();
	private Set<VCFHeaderLine> gvcfHeaderMetaInfo;
	private JointCallingPrepareOptions options = null;
	private GenomeLocationParser parser = null;
	private ReferenceShare genomeShare = null;
	private DbsnpShare dbsnpShare = null;
	private VCFLocalLoader loader = null;
	private VariantRegionFilter filter = null;
	private MultipleVCFHeaderForJointCalling headers = new MultipleVCFHeaderForJointCalling();
	private LazyVCFGenotypesContext.HeaderDataCache vcfHeaderDataCache =
    		new LazyVCFGenotypesContext.HeaderDataCache();
	SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
	String formatStr =formatter.format(new Date());
	private Boolean oneByOne=true;
	Set<String> chrs=new LinkedHashSet<>();
	Map<String,Integer> chr_region_start=new LinkedHashMap<>();
	Map<String,Integer> chr_region_end=new LinkedHashMap<>();
	VCFCodec query_codec=new VCFCodec(); 
	VCFHeaderVersion version = null;
	private VariantAnnotatorEngine annotationEngine;
	protected List<String> annotationsToUse = new ArrayList<>();
	
	protected List<String> annotationGroupsToUse = new ArrayList<>(
			Arrays.asList(new String[] { StandardAnnotation.class.getSimpleName() }));
	private ArrayList<String> mapGvcfList=new ArrayList<String>();
	private ArrayList<FeatureReader<VariantContext> > gvcf_reader_array=new ArrayList<FeatureReader<VariantContext> >();
	private ArrayList<BufferedReader> gvcfBufReader=new ArrayList<BufferedReader>();
	private ArrayList<VariantContext> curSamplesVC=new ArrayList<VariantContext>();
	private ArrayList<String> samplesName=new  ArrayList<String>();
	private HashMap<Integer,Set<String>> mapIndexSample=new HashMap<>();
	Map<String,String> pathSample=new HashMap();
	HashMap<String, Integer> sampleIndex=new HashMap();
	Comparator<VariantContext> comparator4 = new Comparator<VariantContext>() {
		public int compare(VariantContext t1,VariantContext t2) {
			return t1.getStart()-t2.getStart();
		}
	};
//	CloseableTribbleIterator<VariantContext>[] vc_iter_array=new CloseableTribbleIterator[patch_samples];
	@Override
	protected void setup(Context context) throws IOException, InterruptedException {//setup在map函数前运行，只运行一次
		conf = context.getConfiguration();
		options = new JointCallingPrepareOptions();
		options.getOptionsFromHadoopConf(conf);
		System.out.println(accumulate.size());
		FileSplit fileSplit = (FileSplit)context.getInputSplit();
		final Path list_file=fileSplit.getPath();
		final FileSystem fs = list_file.getFileSystem(context.getConfiguration());
		final FSDataInputStream ins = fs.open(list_file);
		AsciiLineReader splitReader = new AsciiLineReader(ins);
		AsciiLineReaderIterator split_it = new AsciiLineReaderIterator(splitReader);
		long split_end=fileSplit.getStart()+fileSplit.getLength()-1;
		tmpOutStr=String.valueOf(fileSplit.getStart());
		if(fileSplit.getStart()==0) {
			ins.seek(0);
		}else {
			ins.seek(fileSplit.getStart());
			ins.readLine();
		}
		System.out.println("GVCF files assigned to this mapper:");
		while(ins.getPos()<=split_end) {
			String tmp_line=ins.readLine();
			if(tmp_line==null) {
				break;
			}
			mapGvcfList.add(tmp_line);
			System.out.println("\t"+tmp_line);
		}
		processed_samples=mapGvcfList.size();
		split_it.close();
		splitReader.close();
		ins.close();
		Set<String> mapSamplePath=new HashSet<>();
		for(String gvcfPath:mapGvcfList) {
			Path path2=new Path(gvcfPath);
			String[] eles=gvcfPath.split("/");
			mapSamplePath.add(eles[eles.length-1]);
		}
		
		String sampleIndexFile=options.getOutput()+"/../vcfheaderinfo";
		Path sampleIndexFilePath=new Path(sampleIndexFile);
		BufferedReader indexReader = new BufferedReader(new InputStreamReader(sampleIndexFilePath.getFileSystem(conf).open(sampleIndexFilePath)));
		String indexLine;
		Logger logger=LoggerFactory.getLogger(JointCallingPrepareMapper.class);
		
		logger.warn("before hash");
		while((indexLine=indexReader.readLine())!=null) {
			String[] eles=indexLine.split("\t");
			if(eles.length!=3) {
				logger.error("vcfheaderinfo file format error");
			}
			String name;
			if(eles[2].endsWith(",")) {
				name=eles[2].substring(0,eles[2].length()-1);
			}else {
				name=eles[2];
			}
			if(mapSamplePath.contains(eles[0])) {
				sampleIndex.put(name, Integer.parseInt(eles[1]));
				pathSample.put(eles[0], name);
			}
		}
		mapSamplePath.clear();
		logger.warn("pathSample size:\t"+pathSample.size());
		indexReader.close();


		realMapperNumber=Integer.parseInt(conf.get(JointCallingPrepare.Real_Mappers));
		Path path = new Path(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP));
		//SeekableStream in = WrapSeekable.openPath(vVCF.getFileSystem(conf), vVCF);
		SeekableStream in=new SeekableFileStream(new File(options.getTmpOut()+"/virtual.vcf"));
		merged_header = VCFHeaderReader.readHeaderFrom(in);//从header文件中读取header
		logger.warn("after get merged header");
		in.close();
		gvcfHeaderMetaInfo=merged_header.getMetaDataInInputOrder();
		//Set<VCFHeader> mapMultiSamplesHeaderSet=new HashSet<VCFHeader>();
		Integer sISize=sampleIndex.size();
		for(String gvcfPath:mapGvcfList) {
			Path path2=new Path(gvcfPath);
//			SeekableStream in2 = WrapSeekable.openPath(path2.getFileSystem(conf), path2);
//			VCFHeader tmp_header=VCFHeaderReader.readHeaderFrom(in2);
			String[] eles=gvcfPath.split("/");
			if(!pathSample.containsKey(eles[eles.length-1])) {
				logger.error("no such path in vcfHeaderInfo");
			}
			String sampleName=pathSample.get(eles[eles.length-1]);
			sampleNames.add(sampleName);
			mapSMtagInt+=sampleIndex.get(sampleName);
			if(version==null) {
				for (final VCFHeaderLine line : gvcfHeaderMetaInfo) {
					if (VCFHeaderVersion.isFormatString(line.getKey())) {
						version = VCFHeaderVersion.toHeaderVersion(line.getValue());
						break;
					}
				}
			}
			
			//VCFCodec tmp_codec=new VCFCodec();
			Set<String> curSample=new HashSet<>();
			curSample.add(sampleName);
			mapIndexSample.put(sampleIndex.get(sampleName),curSample);
			//tmp_codec.setVCFHeader(new VCFHeader(gvcfHeaderMetaInfo,curSample), version);
			HeaderDataCache vcfHeaderDataCache = new HeaderDataCache();
//			codec.add(tmp_codec);
			if(gvcfPath.startsWith("file://")) {
				gvcfPath=gvcfPath.substring(7);
			}
		}
		logger.warn("sampleNames Size:\t"+sampleNames.size());
		logger.warn("mapGvcfList Size:\t"+mapGvcfList.size());
		mapMergedHeader=new VCFHeader(gvcfHeaderMetaInfo,sampleNames);
		//mapMergedHeader=getVCFHeaderFromInput(mapMultiSamplesHeaderSet);
		
		if(merged_header == null)
			throw new RuntimeException("header is null !!!");
		logger.warn("after get merged header");
		contigs = new HashMap<>();
		chrIndexs = new HashMap<>();
		for (VCFContigHeaderLine line : merged_header.getContigLines()) {
			chrIndexs.put(line.getID(), line.getContigIndex());
		}
		for (VCFContigHeaderLine line : merged_header.getContigLines()) {
			contigs.put(line.getContigIndex(), line.getID());
		}
		windowSize = options.getWindowsSize();


		for (Map.Entry<Integer, String> entry : contigs.entrySet()) {
			// System.out.println(entry.getKey()+"\t"+entry.getValue());
			String chr = entry.getValue();
			int contigLength = merged_header.getSequenceDictionary().getSequence(chr).getSequenceLength();
			accumulate.put(entry.getKey(), refLength);
			refLength += contigLength;
		}
	}
	
	@Override
	protected void map(LongWritable key, Text value, Context context)
			throws IOException, InterruptedException {
		//#START each mapper is assigned some gvcf files, do the combineGvcfs, then genotypeGvcfs in reduce step
			// traverse the genome in a certain window size;
		//ArrayList<VariantContext> preUn=new ArrayList<VariantContext>();
			//System.out.println(formatter.format(new Date())+"\tmap start\t"+win_line);
		String[] eles=value.toString().split("/");
		String samplePath=eles[eles.length-1];
		if(!pathSample.containsKey(samplePath)) {
			logger.error("no such path in vcfHeaderInfo");
		}
		
		VCFCodec tmp_codec=new VCFCodec();
		tmp_codec.setVCFHeader(new VCFHeader(gvcfHeaderMetaInfo,mapIndexSample.get(sampleIndex.get(pathSample.get(samplePath)))), version);
		Path samplePathPath=new Path(value.toString());
		BufferedReader reader=new BufferedReader(new InputStreamReader(new GZIPInputStream(samplePathPath.getFileSystem(conf).open(samplePathPath))));
		String line;
		String lastChr="chr0";
		int lastPos=0;
		while((line=reader.readLine())!=null) {
			if(line.startsWith("#")) {
				continue;
			}
			final VariantContext v = tmp_codec.decode(line);
			Integer start_=v.getStart();
//			if(start_>1000000 || chrIndexs.get(v.getContig())>1){
//				break;
//			}
			if(v.getNAlleles()>2) {
				Integer chrIndex=chrIndexs.get(v.getContig());
				String outString=chrIndex+"\t"+v.getStart()+"\t"+v.getEnd();
				Text outvalue=new Text(outString);
				int sWin = v.getStart() / windowSize ;
				int eWin = v.getEnd() / windowSize;
				for (int j = sWin; j <= eWin; j++) {
//					System.out.println(outKey+"\t"+outString);
					outKey.set(chrIndex, j, v.getStart());
					context.write(outKey, outvalue);
				}
				if(!lastChr.equals(v.getContig())){
					System.out.println("cur process chr:\t"+v.getContig());
					lastChr=v.getContig();
					lastPos=0;
				}
				if(lastPos<v.getStart()/10000000){
					System.out.println("cur process pos:\t"+v.getStart());
					lastPos=v.getStart()/10000000;
				}
			}


		}
	}
}