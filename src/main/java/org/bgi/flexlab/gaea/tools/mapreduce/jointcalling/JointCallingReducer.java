package org.bgi.flexlab.gaea.tools.mapreduce.jointcalling;

import java.io.BufferedReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.data.structure.dbsnp.DbsnpShare;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.variant.filter.VariantRegionFilter;
import org.bgi.flexlab.gaea.tools.jointcalling.JointCallingEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

public class JointCallingReducer
		extends Reducer<WindowsBasedWritable, VariantContextWritable, NullWritable, VariantContextWritable> {

	private int windowSize;
	private HashMap<Integer, String> contigs = null;
	private JointCallingEngine engine = null;
	private VariantContextWritable outValue = new VariantContextWritable();

	private JointCallingOptions options = null;
	private GenomeLocationParser parser = null;
	private ReferenceShare genomeShare = null;
	private DbsnpShare dbsnpShare = null;
	private VCFLocalLoader loader = null;
	private VariantRegionFilter filter = null;
	private VCFHeader header = null;
	private VCFHeader header2=null;
	private MultipleVCFHeaderForJointCalling headers = new MultipleVCFHeaderForJointCalling();
	private ArrayList<ArrayList<String> > multiMapSampleNames=new ArrayList<ArrayList<String> >();
	public BufferedReader bp_reader=null;
	public String winLine=null;
	SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
	@Override
	protected void setup(Context context) throws IOException {
		Configuration conf = context.getConfiguration();
		contigs = new HashMap<>();
		options = new JointCallingOptions();
		options.getOptionsFromHadoopConf(conf);
		Path path = new Path(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP));
		SeekableStream in = WrapSeekable.openPath(path.getFileSystem(conf), path);
//		SeekableStream in=new SeekableFileStream(new File(options.getTmpOut()+"/virtual.vcf"));
		System.out.println(formatter.format(new Date())+"\tbefore readHeader");
		header = VCFHeaderReader.readHeaderFrom(in);
		System.out.println(formatter.format(new Date())+"\tafter readHeader");
		in.close();
		if(header == null)
			throw new RuntimeException("header is null !!!");
		
		for (VCFContigHeaderLine line : header.getContigLines()) {
			contigs.put(line.getContigIndex(), line.getID());
		}
		windowSize = options.getWindowsSize();
		parser = new GenomeLocationParser(header.getSequenceDictionary());
		headers.readHeaders(conf);
		String sampleStr = conf.get(JointCalling.INPUT_ORDER);;
		String[] allSample=sampleStr.split(",");
		int mapperLine=allSample.length/options.getMapperNumber()<2?2:allSample.length/options.getMapperNumber();
		ArrayList<String> ele=new ArrayList<String>();
		
		for(int i=0;i<allSample.length;i++) {
			ele.add(allSample[i]);
			if((i+1)%mapperLine==0) {
				ArrayList<String> tmp_ele=new ArrayList<String>();
				tmp_ele.addAll(ele);
				Collections.sort(tmp_ele);
				multiMapSampleNames.add(tmp_ele);
				ele.clear();
			}
		}
		if(ele.size()>0) {
			ArrayList<String> tmp_ele=new ArrayList<String>();
			tmp_ele.addAll(ele);
			Collections.sort(tmp_ele);
			multiMapSampleNames.add(tmp_ele);
		}
		System.out.println(formatter.format(new Date())+"\tbefore engine init");
//		System.out.println("sampleStr:\t"+sampleStr);


		if(sampleStr != null) {
			if(options.getMapperMode()) {
				//JointCallingPrepareOptions options, GenomeLocationParser parser, VCFHeader vcfheader,
				//MultipleVCFHeaderForJointCalling multiHeaders,String[] sampleArray,ArrayList<ArrayList<String> >multiMapSamples,Configuration conf
				engine = new JointCallingEngine(options, parser,header,headers, multiMapSampleNames,conf);
			}else {
				engine = new JointCallingEngine(options, parser,headers);
			}
		}else {
			engine = new JointCallingEngine(options, parser,headers);
		}
		System.out.println(formatter.format(new Date())+"\tafter engine init");
		genomeShare = new ReferenceShare();
		genomeShare.loadChromosomeList(options.getReference());
		dbsnpShare = new DbsnpShare(options.getDBSnp(), options.getReference());
		dbsnpShare.loadChromosomeList(options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
		loader = new VCFLocalLoader(options.getDBSnp());
		filter = new VariantRegionFilter();
		header2=engine.getVCFHeader();
		if(header2 == null)
			throw new RuntimeException("header is null !!!");
		String bpFile=options.getTmpOut()+"/AllBPs";
//		bp_reader=new BufferedReader(new FileReader(bpFile));
//		winLine=bp_reader.readLine();
//		System.out.println("reduce setup done");
	}

	@Override
	public void reduce(WindowsBasedWritable key, Iterable<VariantContextWritable> values, Context context)
			throws IOException, InterruptedException {
		System.out.println(formatter.format(new Date())+"\treduce start");
		int winNum = key.getWindowsNumber();
		int start = winNum * windowSize;
		if(start == 0)
			start = 1;
		String chr = contigs.get(key.getChromosomeIndex());
		int chrInx=key.getChromosomeIndex();
		int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
		int end = Math.min(contigLength, start + windowSize - 1);
		
		long startPosition = dbsnpShare.getStartPosition(chr, winNum, options.getWindowsSize());
		ArrayList<VariantContext> dbsnps = null;
		if(startPosition >= 0)
			dbsnps = filter.loadFilter(loader, chr, startPosition, end);
		engine.init(dbsnps);
		Set<Integer> bps=new TreeSet();


		//原方式
//		int iterTmp=0;
//		for(int iter=start;iter<=end;iter++){
//			VariantContext variantContext = engine.variantCalling2(values.iterator(),
//					parser.createGenomeLocation(chr, iter), genomeShare.getChromosomeInfo(chr));
//			if (variantContext == null){
//				continue;
//			}
//			if(!variantContext.getContig().equals("chr1")){
//				iterTmp++;
//				if(iterTmp<3){
//					System.out.println(variantContext.getContig()+"\t"+variantContext.getStart());
//				}
//			}
//			CommonInfo info = variantContext.getCommonInfo();
//			HashMap<String, Object> maps = new HashMap<>();
//			maps.putAll(info.getAttributes());
//			maps.remove("SM");
//			info.setAttributes(maps);
//			outValue.set(variantContext, header2);
//			context.write(NullWritable.get(), outValue);
//		}

		//为了避免每次从头读取winBp，保留bp_reader，使其只读一遍
		System.out.println(formatter.format(new Date())+"\tbp window before get bps:\t"+winLine);
		while(true) {
			String[] bpeles=winLine.split("\t");
			Integer bpcontig=Integer.valueOf(bpeles[0]);
			Integer bpstart=Integer.valueOf(bpeles[1]);
			Integer bpend=Integer.valueOf(bpeles[2]);
			if(bpcontig==chrInx) {
				if(bpstart<=end && bpend>=start) {
					Integer realStart=bpstart>start?bpstart:start;
					Integer realEnd=bpend<end?bpend:end;
					for(int i=realStart;i<=realEnd;i++) {
						bps.add(i);
					}
					if(bpend<end) {
						if ((winLine = bp_reader.readLine()) == null) {
							return;
						}
					}else{
						break;
					}
				}else if(bpend<start) {
					if((winLine=bp_reader.readLine())==null){
						return;
					}
				}else if(bpstart>end){
					break;
				}
			}else if(bpcontig<chrInx) {
				if((winLine=bp_reader.readLine())==null){
					return;
				}
			}else {
				break;
			}
		}
		System.out.println(formatter.format(new Date())+"\tcurrent reduce key:\t"+chr+"\t"+chrInx+"\t"+start+"\t"+end);
		System.out.println(formatter.format(new Date())+"\tbp window after get bps:\t"+winLine);
		System.out.println(formatter.format(new Date())+"\tbps size in region:\t"+bps.size());
//		for (int iter = start; iter <= end; iter++){
		int iterTmp=0;
		for(Integer iter:bps) {
			VariantContext variantContext = engine.variantCallingReduce(values.iterator(),
					parser.createGenomeLocation(chr, iter), genomeShare.getChromosomeInfo(chr));
			if (variantContext == null){
				continue;
			}
//			iterTmp++;
//			if(iterTmp<3){
//				System.out.println(variantContext);
//				System.out.println(variantContext.getAttribute("BaseQRankSum"));
//			}
			if(!variantContext.getContig().equals("chr1")){
				iterTmp++;
				if(iterTmp<3){
					System.out.println(variantContext.getContig()+"\t"+variantContext.getStart());
				}
			}
			CommonInfo info = variantContext.getCommonInfo();
			HashMap<String, Object> maps = new HashMap<>();
			maps.putAll(info.getAttributes());
			maps.remove("SM");
			info.setAttributes(maps);
			outValue.set(variantContext, header2);
			context.write(NullWritable.get(), outValue);
		}


		//breakPoints.clear();
		//bp_reader.close();
	}


	@Override
	protected void cleanup(Context context) throws IOException{
		bp_reader.close();
	}
}
