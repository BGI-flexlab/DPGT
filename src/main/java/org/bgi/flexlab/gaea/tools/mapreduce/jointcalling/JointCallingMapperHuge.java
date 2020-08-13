package org.bgi.flexlab.gaea.tools.mapreduce.jointcalling;

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

public class JointCallingMapperHuge extends
	Mapper<LongWritable, Text, WindowsBasedWritable, VariantContextWritable>{
	
	//private int windowSize = 10000;
	//private int thread_num=10;
	private Logger logger = LoggerFactory.getLogger(JointCallingMapperHuge.class);
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
	private JointCallingOptions options = null;
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
	Comparator<VariantContext> comparator4 = new Comparator<VariantContext>() {
		public int compare(VariantContext t1,VariantContext t2) {
			return t1.getStart()-t2.getStart();
		}
	};
	public static HashMap<Integer,Long> accumulate=new HashMap<>();
//	CloseableTribbleIterator<VariantContext>[] vc_iter_array=new CloseableTribbleIterator[patch_samples];
	@Override
	protected void setup(Context context) throws IOException, InterruptedException {//setup在map函数前运行，只运行一次
		conf = context.getConfiguration();
		options = new JointCallingOptions();
		options.getOptionsFromHadoopConf(conf);
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
		HashMap<String, Integer> sampleIndex=new HashMap();
		String sampleIndexFile=options.getOutput()+"/../vcfheaderinfo";
		Path sampleIndexFilePath=new Path(sampleIndexFile);
		BufferedReader indexReader = new BufferedReader(new InputStreamReader(sampleIndexFilePath.getFileSystem(conf).open(sampleIndexFilePath)));
		String indexLine;
		Logger logger=LoggerFactory.getLogger(JointCallingMapperHuge.class);
		logger.warn("in mapperHuge");
		Map<String,String> pathSample=new HashMap();
		int totalSampleSize=0;
		while((indexLine=indexReader.readLine())!=null) {
			totalSampleSize++;
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
		indexReader.close();
		realMapperNumber=Integer.parseInt(conf.get(JointCalling.Real_Mappers));
		//Path path = new Path(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP));
		//SeekableStream in = WrapSeekable.openPath(path.getFileSystem(conf), path);
		SeekableStream in=new SeekableFileStream(new File(options.getTmpOut()+"/virtual.vcf"));
		merged_header = VCFHeaderReader.readHeaderFrom(in);//从header文件中读取header
		in.close();
		gvcfHeaderMetaInfo=merged_header.getMetaDataInInputOrder();
		//Set<VCFHeader> mapMultiSamplesHeaderSet=new HashSet<VCFHeader>();
		Integer sISize=sampleIndex.size();
		int i=0;
		int mapGvcfListIndex=0;
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
			if(mapGvcfListIndex==0) {
				mapSMtagInt = sampleIndex.get(sampleName);
			}
			mapGvcfListIndex++;
			if(version==null) {
				for (final VCFHeaderLine line : gvcfHeaderMetaInfo) {
					if (VCFHeaderVersion.isFormatString(line.getKey())) {
						version = VCFHeaderVersion.toHeaderVersion(line.getValue());
						break;
					}
				}
			}
			
			VCFCodec tmp_codec=new VCFCodec();
			Set<String> curSample=new HashSet<>();
			curSample.add(sampleName);
			mapIndexSample.put(i,curSample);
			tmp_codec.setVCFHeader(new VCFHeader(gvcfHeaderMetaInfo,curSample), version);
			HeaderDataCache vcfHeaderDataCache = new HeaderDataCache();
			codec.add(tmp_codec);
			logger.warn(i+"\tcurrent index");
			if(gvcfPath.startsWith("file://")) {
				gvcfPath=gvcfPath.substring(7);
			}
			BufferedReader reader=null;
			if(gvcfPath.endsWith(".gz")) {
				if(gvcfPath.startsWith("/user")) {
					reader=new BufferedReader(new InputStreamReader(new GZIPInputStream(path2.getFileSystem(conf).open(path2))));
				}else {
					reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gvcfPath))));
				}
			}else {
				reader = new BufferedReader(new InputStreamReader(new FileInputStream(gvcfPath)));
			}
			gvcfBufReader.add(reader);
			i++;
		}
		logger.warn("sampleNames Size:\t"+sampleNames.size());
		logger.warn("mapGvcfList Size:\t"+mapGvcfList.size());
		pathSample.clear();
		sampleIndex.clear();
		logger.warn("totalSampleSize:\t"+totalSampleSize+"\tmapSMtagInt:\t"+mapSMtagInt);
		mapSMtagInt+=totalSampleSize;
		System.out.println("mapSMtag:\t"+mapSMtagInt);
		mapMergedHeader=new VCFHeader(gvcfHeaderMetaInfo,sampleNames);
		//mapMergedHeader=getVCFHeaderFromInput(mapMultiSamplesHeaderSet);
		vcfHeaderDataCache.setHeader(mapMergedHeader);
		
		mapSMtag=Utils.join("_",mapMergedHeader.getSampleNamesInOrder());
		if(mapMergedHeader.getSampleNamesInOrder().size()==1) {
			mapSMtag=mapSMtag+"_";
		}
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
		Long refLength=0L;
		for (Map.Entry<Integer, String> entry : contigs.entrySet()) {
			// System.out.println(entry.getKey()+"\t"+entry.getValue());
			String chr = entry.getValue();
			int contigLength = merged_header.getSequenceDictionary().getSequence(chr).getSequenceLength();
			accumulate.put(entry.getKey(), refLength);
			refLength += contigLength;
		}
		WindowsBasedWritable winBW=new WindowsBasedWritable(accumulate,0L,refLength,options.getReducerNumber());
		windowSize = options.getWindowsSize();
		parser = new GenomeLocationParser(merged_header.getSequenceDictionary());
		
		logger.warn("after parser");
		
		String sampleStr = conf.get(JointCalling.INPUT_ORDER,null);
//		if(sampleStr != null)
//			engine = new JointCallingEngine(options, parser,headers);
//		else
//			engine = new JointCallingEngine(options, parser,headers);
		genomeShare = new ReferenceShare();
		genomeShare.loadChromosomeList(options.getReference());
		dbsnpShare = new DbsnpShare(options.getDBSnp(), options.getReference());
		dbsnpShare.loadChromosomeList(options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
		loader = new VCFLocalLoader(options.getDBSnp());
		filter = new VariantRegionFilter();
		for (final VCFContigHeaderLine contig : merged_header.getContigLines())
			contigDict.put(contig.getID(), contig.getContigIndex());
		logger.warn("after contigDict Map done");
		//System.out.println(merged_header);
		//vcfHeaderDataCache.setHeader(merged_header);
		annotationEngine = new VariantAnnotatorEngine(annotationGroupsToUse, annotationsToUse,
				Collections.<String>emptyList());
		annotationEngine.initializeDBs(options.getDBSnp() != null);
		logger.warn("setup done");
		// get merged region from win split
		
	}
	
	@Override
	protected void map(LongWritable key, Text value, Context context)
			throws IOException, InterruptedException {
		//#START each mapper is assigned some gvcf files, do the combineGvcfs, then genotypeGvcfs in reduce step
			// traverse the genome in a certain window size;
		run_times++;
		if(run_times>1) {
			return;
		}

		Integer number=0;
		Integer lastchrIndex=-1;
		//BufferedReader win_reader=new BufferedReader(new FileReader(conf.get(JointCallingPrepare.Window_File)));
		String winFilePath=conf.get(JointCalling.Window_File);
		//FileSystem fs=FileSystem.get(URI.create(winFilePath),conf);
		if(winFilePath.startsWith("file://")) {
			winFilePath=winFilePath.substring(7);
		}
		BufferedReader win_reader=new BufferedReader(new FileReader(winFilePath));
		String win_line=null;
		gvcfBufReader.clear();
		Path firstPath=new Path(mapGvcfList.get(0));
		FileSystem gvcf_fs=firstPath.getFileSystem(conf);
		for(int i=0;i!=mapGvcfList.size();i++) {
			String gvcfFile=mapGvcfList.get(i);
			if(gvcfFile.startsWith("file://")) {
				gvcfFile=gvcfFile.substring(7);
			}
			BufferedReader reader=null;
			if(gvcfFile.endsWith(".gz")) {
				if(gvcfFile.startsWith("/user")) {
					reader=new BufferedReader(new InputStreamReader(new GZIPInputStream(gvcf_fs.open(new Path(gvcfFile)))));
				}else {
					reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gvcfFile))));
				}
			}else {
				reader = new BufferedReader(new InputStreamReader(new FileInputStream(gvcfFile)));
			}
			gvcfBufReader.add(reader);
		}
		String[] tmpLines=new String[mapGvcfList.size()];
		for(int i=0;i!=mapGvcfList.size();i++) {
			String tmpline=gvcfBufReader.get(i).readLine();
			tmpLines[i]=tmpline;
		}
		//method 3, merge code from MapperHuge
		win_reader=new BufferedReader(new FileReader(winFilePath));
		win_line=null;
		String last_contig=null;
		ChromosomeInformationShare ref=null;
		int vNum=0;
		String bpFile=options.getTmpOut()+"/AllBPs";
		BufferedReader bp_reader=new BufferedReader(new FileReader(bpFile));
		String winLine=bp_reader.readLine();
		String[] bpeles=winLine.split("\t");
		Integer bpcontig=Integer.valueOf(bpeles[0]);
		Integer bpstart=Integer.valueOf(bpeles[1]);
		Integer bpend=Integer.valueOf(bpeles[2]);
		ArrayList<Integer> bpStarts=new ArrayList<Integer>();
		ArrayList<Integer> bpEnds=new ArrayList<Integer>();
		Boolean bpFileEnd=false;
		Boolean logOut=false;
		Integer lastLogPoint=0;
		Integer last_end=0;
		Integer last_pos=0;
		Integer last_real_pos=-1;
		Integer last_bp_start=-1;
		Integer last_bp_end=-1;
		Set<Integer> last_end_breakpoints=new TreeSet<Integer>();
		while((win_line=win_reader.readLine())!=null) {
			if(logOut)
				System.out.println(formatter.format(new Date())+"\tmap start\t"+win_line);
			String[] eles=win_line.split("\t");
			if(eles.length!=3) {
				System.err.println("the format of input windows file errors!");
				System.exit(-1);
			}
			String contig=eles[0];
			Integer start=Integer.valueOf(eles[1]);
			Integer end=Integer.valueOf(eles[2]);
			if(start>10000000){
				break;
			}
			String chr = contig;
			Integer curLogPoint=(int)start/5000000;
			if(!curLogPoint.equals(lastLogPoint)) {
				logOut=true;
			}else {
				logOut=false;
			}
			lastLogPoint=curLogPoint;
			Integer winChrID=chrIndexs.get(contig);
			if(last_bp_start!=-1){
				bpStarts.add(start);
				bpEnds.add(last_bp_end);
				last_bp_start=-1;
				last_bp_end=-1;
			}
			while(!bpFileEnd) {
				if(bpcontig.equals(winChrID)) {
					if(bpstart<=end && bpend>=start) {
						bpStarts.add(bpstart);
						bpEnds.add(bpend);
						winLine=bp_reader.readLine();
						if(winLine==null) {
							bpFileEnd=true;
							bp_reader.close();
							break;
						}
						bpeles=winLine.split("\t");
						bpcontig=Integer.valueOf(bpeles[0]);
						bpstart=Integer.valueOf(bpeles[1]);
						bpend=Integer.valueOf(bpeles[2]);
					}else if(bpstart>end) {
						if(bpend>end){
							last_bp_start=bpstart;
							last_bp_end=bpend;
						}
						break;
					}else {
						winLine=bp_reader.readLine();
						if(winLine==null) {
							bpFileEnd=true;
							bp_reader.close();
							break;
						}
						bpeles=winLine.split("\t");
						bpcontig=Integer.valueOf(bpeles[0]);
						bpstart=Integer.valueOf(bpeles[1]);
						bpend=Integer.valueOf(bpeles[2]);
					}
				}else if(bpcontig>winChrID) {
					break;
				}else {
					winLine=bp_reader.readLine();
					if(winLine==null) {
						bpFileEnd=true;
						bp_reader.close();
						break;
					}
					bpeles=winLine.split("\t");
					bpcontig=Integer.valueOf(bpeles[0]);
					bpstart=Integer.valueOf(bpeles[1]);
					bpend=Integer.valueOf(bpeles[2]);
				}
			}
			LinkedList<VariantContext> region_vcs=new LinkedList<VariantContext>();
			Set<Integer> end_breakpoints=new TreeSet<Integer>();
			if(logOut) {
				System.out.println(formatter.format(new Date())+"\tbefore query");
			}
			for(int i=0;i!=curSamplesVC.size();i++) {
				VariantContext tmpVC=curSamplesVC.get(i);
				if(tmpVC.getChr().equals(chr) && tmpVC.getStart()<=end && tmpVC.getEnd()>=start) {
					region_vcs.add(tmpVC);
					if(tmpVC.getNAlleles()>2) {
						for(int tmpPos=tmpVC.getStart();tmpPos<=tmpVC.getEnd();tmpPos++) {
							end_breakpoints.add(tmpPos);
						}
					}else {
						end_breakpoints.add(tmpVC.getEnd());
					}
					if(tmpVC.getEnd()<=end) {
						curSamplesVC.remove(i);
						i--;
					}
				}
			}
			for(int i=0;i!=mapGvcfList.size();i++) {
				String line=null;
				String sampleName=sampleNames.get(i);
				int index=0;
				while((line=gvcfBufReader.get(i).readLine())!=null) {
					if(line.startsWith("#")) {
						continue;
					}
					final VariantContext v = codec.get(i).decode(line);
					
					CommonInfo info = v.getCommonInfo();
					if(!info.hasAttribute("SM"))
						info.putAttribute("SM",sampleName);
					Boolean include=false;
					if(v.getChr().equals(chr)) {
						if(v.getStart()<=end && v.getEnd()>=start) {
							if(v.getNAlleles()>2) {
								for(int tmpPos=v.getStart();tmpPos<=v.getEnd();tmpPos++) {
									end_breakpoints.add(tmpPos);
								}
							}else {
								end_breakpoints.add(v.getEnd());
							}
							region_vcs.add(v);
							if(v.getEnd()>end) {
								curSamplesVC.add(v);
							}
						}else if(v.getStart()>end) {
							curSamplesVC.add(v);
							break;
						}else if(v.getEnd()<start) {
							continue;
						}else {
						}
					}else if(contigDict.get(v.getChr()) > contigDict.get(chr) ){
						curSamplesVC.add(v);
						break;
					}else{
						continue;
					}
				}
			}
			if(logOut) {
				System.out.println(formatter.format(new Date())+"\tafter query");
			}
			Collections.sort(region_vcs,comparator4);
			//do combineGVCFS
			
			if(!contig.equals(last_contig)) {
				ref=genomeShare.getChromosomeInfo(contig);
				last_pos=0;
			}else {
				last_pos=start;
			}
			ArrayList<Integer> realBreakpointStart=new ArrayList<Integer>();
			ArrayList<Integer> realBreakpointEnd=new ArrayList<Integer>();
			int bpIndex=0;
			for(Integer pos:end_breakpoints) {
				Boolean whether_keep=false;
				for(;bpIndex!=bpStarts.size();bpIndex++) {
					if(last_pos<=bpEnds.get(bpIndex) && pos>=bpStarts.get(bpIndex)) {
						whether_keep=true;
						break;
					}else if(last_pos>bpEnds.get(bpIndex)) {
						continue;
					}else if(pos<bpStarts.get(bpIndex)) {
						break;
					}
				}
				if(whether_keep) {
					realBreakpointStart.add(last_pos);
					realBreakpointEnd.add(pos);
					if(pos>end){
						last_end_breakpoints.add(pos);
					}
					last_real_pos=pos+1;
				}
				last_pos=pos+1;
			}

			if(region_vcs.size()>0) {
				if(logOut) {
					System.out.println(formatter.format(new Date())+"\tvc number in this region:\t"+region_vcs.size()+"\tcurSamplesVC:\t"+curSamplesVC.size()+"\tbpStarts:\t"+bpStarts.size()+"\tend_breakpoints:\t"+end_breakpoints.size()+"\trealBreakpoints:\t"+realBreakpointStart.size());
				}
			}
			end_breakpoints.clear();
			bpStarts.clear();
			bpEnds.clear();
			if(!contig.equals(last_contig)) {
				last_end=0;
			}
			//last_pos=0;
			int tmpIter=0;
			for(int i=0;i!=realBreakpointStart.size();i++) {
				Integer posStart=realBreakpointStart.get(i);
				Integer pos=realBreakpointEnd.get(i);

				if(pos>end) {
					break;
				}
				if(pos<start) {
					continue;
				}
				List<VariantContext> stoppedVCs = new ArrayList<>();
				//System.out.println(region_vcs.size());
				for(int j=0;j<region_vcs.size();j++) {
					final VariantContext vc = region_vcs.get(j);
					//the VC for the previous state will be stopped if its position is previous to the current position or it we've moved to a new contig
		            if ( vc.getStart() <= pos) {
		            	
		            	if(vc.getEnd()<posStart) {
		            		region_vcs.remove(j);
		            		j--;
		            	}else {
		            		stoppedVCs.add(vc);
		            	}
		                // if it was ending anyways, then remove it from the future state
		                if ( vc.getEnd() == pos) {
		                	//region_vcws.samples.removeAll(vc.getSampleNames());
		                	if(j>=0) {
		                		region_vcs.remove(j);
		                		j--;
		                	}
		                }
		            }else {
		            	break;
		            }
				}
				//System.out.println(region_vcs.size());
				if ( !stoppedVCs.isEmpty()) {
					//System.out.println(stoppedVCs.size());
//					for(VariantContext vc:stoppedVCs) {
//						System.out.println(vc);
//					}
		            final GenomeLocation gLoc =parser.createGenomeLocation(chr, last_end);
		            final VariantContext mergedVC;
		            final Byte refBase =(byte) ref.getBase(gLoc.getStart());
		            int curStart=last_end+1;
		            if ( containsTrueAltAllele(stoppedVCs) ) {
		            	vNum++;
		                mergedVC = ReferenceConfidenceVariantContextMerger.mapMerge(stoppedVCs, parser.createGenomeLocation(chr, posStart), (byte) ref.getBase(posStart-1), false, false, annotationEngine);
						tmpIter++;
//						if(mergedVC.getStart()==1002417){
//							System.out.println(mergedVC+"\n"+mergedVC.getAttribute("DP"));
//						}
		            }else {
		                mergedVC = referenceBlockMerge(stoppedVCs, gLoc,refBase, pos);
//						if(mergedVC.getStart()==1002417){
//							System.out.println(mergedVC+"\n"+mergedVC.getAttribute("DP"));
//						}
		            }
		            if(mergedVC==null) {
		            	continue;
					}

		            CommonInfo info = mergedVC.getCommonInfo();
		    		if(!info.hasAttribute("SM"))
		    			info.putAttribute("SM",mapSMtagInt);
		    		HashMap<String, Object> maps = new HashMap<>();
		            maps.putAll(info.getAttributes());
		    		info.setAttributes(maps);
					VariantContextWritable outvalue=new VariantContextWritable();
					outvalue.set(mergedVC, mapMergedHeader);
					int sWin = mergedVC.getStart() / windowSize ;
					int eWin = mergedVC.getEnd() / windowSize;
					int chrIndex = chrIndexs.get(mergedVC.getContig());
					if(lastchrIndex==chrIndex){
						number++;
					}else{
						System.out.println(lastchrIndex+"\t"+number);
						number=1;
					}
					lastchrIndex=chrIndex;
					for (int j = sWin; j <= eWin; j++) {
						outKey.set(chrIndex, j, mergedVC.getStart());
						context.write(outKey, outvalue);
					}
					stoppedVCs.clear();
		        }
				last_end=pos;
			}
			last_contig=contig;
			end_breakpoints.clear();
		}
		win_reader.close();
		
		return;
	}
	private boolean containsTrueAltAllele(final List<VariantContext> VCs) {
        if ( VCs == null ) throw new IllegalArgumentException("The list of VariantContexts cannot be null");

        for ( final VariantContext vc : VCs ) {
            if ( vc.getNAlleles() > 2 )
                return true;
        }
        return false;
    }
	private VariantContext referenceBlockMerge(final List<VariantContext> VCs, final GenomeLocation loc,byte refAfter, final int end) {

        final VariantContext first = VCs.get(0);
        
        // ref allele and start
        final Allele refAllele;
        final int start;
        if ( loc == null || !loc.getContig().equals(first.getChr()) || first.getStart() >= loc.getStart() + 1) {
            start = first.getStart();
            refAllele = first.getReference();
        } else {
            start = loc.getStart() + 1;
            refAllele = Allele.create(refAfter, true);
        }
        // attributes
        final Map<String, Object> attrs = new HashMap<>(1);

        // genotypes
        final GenotypesContext genotypes = GenotypesContext.create();
        for ( final VariantContext vc : VCs ) {
            for ( final Genotype g : vc.getGenotypes() )
                genotypes.add(new GenotypeBuilder(g).alleles(GaeaGvcfVariantContextUtils.noCallAlleles(g.getPloidy())).make());
        }
        return new VariantContextBuilder("", first.getChr(), start, end, Arrays.asList(refAllele, VariantContextMerger.NON_REF_SYMBOLIC_ALLELE)).attributes(attrs).genotypes(genotypes).make();
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
	private Set<String> getSampleList(Set<VCFHeader> headers){
		Set<String> samples = new TreeSet<String>();
		for(VCFHeader header : headers){
			for ( String sample : header.getGenotypeSamples() ) {
				samples.add(GaeaGvcfVariantContextUtils.mergedSampleName(null, sample, false));
			}
		}
		return samples;
	}
}