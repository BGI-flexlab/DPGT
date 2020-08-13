package org.bgi.flexlab.gaea.tools.mapreduce.jointcallingprepare;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
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
import org.bgi.flexlab.gaea.util.Utils;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

public class JointCallingPrepareReducer
		extends Reducer<WindowsBasedWritable, Text, NullWritable, Text> {

	private int windowSize;
	private HashMap<Integer, String> contigs = null;
	private JointCallingEngine engine = null;

	private JointCallingPrepareOptions options = null;
	private GenomeLocationParser parser = null;
	private ReferenceShare genomeShare = null;
	private DbsnpShare dbsnpShare = null;
	private VCFLocalLoader loader = null;
	private VariantRegionFilter filter = null;
	private VCFHeader header = null;
	private MultipleVCFHeaderForJointCalling headers = new MultipleVCFHeaderForJointCalling();
	private ArrayList<ArrayList<String> > multiMapSampleNames=new ArrayList<ArrayList<String> >();
	String bpFile=null;
	@Override
	protected void setup(Context context) throws IOException {
		Configuration conf = context.getConfiguration();
		contigs = new HashMap<>();
		
		Path path = new Path(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP));
		SeekableStream in = WrapSeekable.openPath(path.getFileSystem(conf), path);
		header = VCFHeaderReader.readHeaderFrom(in);
		in.close();
		
		if(header == null)
			throw new RuntimeException("header is null !!!");
		
		for (VCFContigHeaderLine line : header.getContigLines()) {
			contigs.put(line.getContigIndex(), line.getID());
		}
		
		
		
		options = new JointCallingPrepareOptions();
		options.getOptionsFromHadoopConf(conf);
		windowSize = options.getWindowsSize();
	}

	@Override
	public void reduce(WindowsBasedWritable key, Iterable<Text> values, Context context)
			throws IOException, InterruptedException {
		int winNum = key.getWindowsNumber();
		int start = winNum * windowSize;
		if(start == 0)
			start = 1;
		String chr = contigs.get(key.getChromosomeIndex());
		int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
		int end = Math.min(contigLength, start + windowSize - 1);
//		System.out.println("here\t"+chr+"\t"+winNum);
		Set<Integer> realBreakpoints=new TreeSet<Integer>();
		int ii=0;
		for (int iter = start; iter <= end; iter++){
			Iterator<Text> it=values.iterator();
			while (it.hasNext()) {
				ii++;
				Text outValue = it.next();
				System.out.println(outValue);
				String eles[]=outValue.toString().split("\t");
				if(eles.length!=3) {
					System.out.println("code error");
					System.exit(-1);
				}
				int vStart=Integer.parseInt(eles[1]);
				int vEnd=Integer.parseInt(eles[2]);
				if(vEnd<=end) {
					for(int i=vStart;i<=vEnd;i++) {
						realBreakpoints.add(i);
					}
				}else {
					for(int i=vStart;i<=end;i++) {
						realBreakpoints.add(i);
					}
				}
			}
		}
//		System.out.println("size:\t"+ii+"\tbreakpoints Size:\t"+realBreakpoints.size());
		int wStart=0;
		int wEnd=0;
		int lastPos=0;
		for(Integer pos:realBreakpoints) {
			if(wStart==0) {
				wStart=pos;
				wEnd=pos;
			}else {
				if(pos-lastPos==1) {
					wEnd=pos;
				}else {
					String outLine=key.getChromosomeIndex()+"\t"+wStart+"\t"+wEnd;
					Text outValue=new Text(outLine);
					context.write(NullWritable.get(), outValue);
//					writer2.write(outLine);
					wStart=pos;
					wEnd=pos;
				}
			}
			lastPos=pos;
		}
		if(wStart!=0) {
			String outLine=key.getChromosomeIndex()+"\t"+wStart+"\t"+wEnd;
			Text outValue=new Text(outLine);
			context.write(NullWritable.get(), outValue);
		}
		realBreakpoints.clear();
	}
}
