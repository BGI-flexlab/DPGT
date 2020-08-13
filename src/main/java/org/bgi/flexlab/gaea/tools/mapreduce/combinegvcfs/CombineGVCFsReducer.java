package org.bgi.flexlab.gaea.tools.mapreduce.combinegvcfs;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.Reducer.Context;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.data.structure.dbsnp.DbsnpShare;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.variant.filter.VariantRegionFilter;
import org.bgi.flexlab.gaea.tools.jointcalling.JointCallingEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCalling;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCallingOptions;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.GaeaVariantContextUtils;

import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

public class CombineGVCFsReducer extends Reducer<WindowsBasedWritable, VariantContextWritable, NullWritable, VariantContextWritable>{
	private int windowSize;
	private HashMap<Integer, String> contigs = null;
	private CombineGVCFsEngine engine = null;
	private VariantContextWritable outValue = new VariantContextWritable();

	private JointCallingOptions options = null;
	private GenomeLocationParser parser = null;
	private ReferenceShare genomeShare = null;
	private DbsnpShare dbsnpShare = null;
	private VCFLocalLoader loader = null;
	private VariantRegionFilter filter = null;
	private VCFHeader header = null;
	private MultipleVCFHeaderForJointCalling headers = new MultipleVCFHeaderForJointCalling();

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

		options = new JointCallingOptions();
		options.getOptionsFromHadoopConf(conf);
		
		windowSize = options.getWindowsSize();
		parser = new GenomeLocationParser(header.getSequenceDictionary());
		headers.readHeaders(conf);
		
		String sampleStr = conf.get(JointCalling.INPUT_ORDER,null);
		if(sampleStr != null)
			engine = new CombineGVCFsEngine(options, parser,header,headers,sampleStr.split(","));
		else
			engine = new CombineGVCFsEngine(options, parser,header,headers,null);
		genomeShare = new ReferenceShare();
		genomeShare.loadChromosomeList(options.getReference());
		dbsnpShare = new DbsnpShare(options.getDBSnp(), options.getReference());
		dbsnpShare.loadChromosomeList(options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
		loader = new VCFLocalLoader(options.getDBSnp());
		filter = new VariantRegionFilter();
		header = engine.getVCFHeader();
		
		if(header == null)
			throw new RuntimeException("header is null!!!");
	}

	@Override
	public void reduce(WindowsBasedWritable key, Iterable<VariantContextWritable> values, Context context)
			throws IOException, InterruptedException {
		int winNum = key.getWindowsNumber();
		int start = winNum * windowSize;
		if(start == 0)
			start = 1;
		String chr = contigs.get(key.getChromosomeIndex());

		int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
		int end = Math.min(contigLength, start + windowSize - 1);

		long startPosition = dbsnpShare.getStartPosition(chr, winNum, options.getWindowsSize());
		ArrayList<VariantContext> dbsnps = null;
		if(startPosition >= 0)
			dbsnps = filter.loadFilter(loader, chr, startPosition, end);
		engine.init(dbsnps);
		
		
		Set<Integer> breakpoints=new TreeSet<>();
		List<VariantContext> tmp_vcs=new ArrayList<>();
		VariantContext tmp_vc;
		
		for(VariantContextWritable it : values) {
			
			tmp_vc=it.get();
			tmp_vcs.add(tmp_vc);
			breakpoints.add(tmp_vc.getStart());
			breakpoints.add(tmp_vc.getEnd());
		}
		int last_pos=0;
		System.out.println(tmp_vcs.size());
		//Iterator<VariantContextWritable> test_it=values.iterator();
		System.out.println("breakpoints number:\t"+breakpoints.size());
		for(Integer it:breakpoints) {
			if(it>end) {
				break;
			}
			if(last_pos==0) {
				last_pos=it;
				continue;
			}else {
				
				VariantContext variantContext=engine.variantCalling3(tmp_vcs.iterator(), parser.createGenomeLocation(chr,last_pos,it-1), genomeShare.getChromosomeInfo(chr));
				if (variantContext == null)
					continue;
				CommonInfo info = variantContext.getCommonInfo();
				HashMap<String, Object> maps = new HashMap<>();
				maps.putAll(info.getAttributes());
				maps.remove("SM");
				info.setAttributes(maps);

				outValue.set(variantContext, header);
				context.write(NullWritable.get(), outValue);
				last_pos=it;
			}
			
		}
//		for (int iter = start; iter <= end; iter++) {
////			if ( containsTrueAltAllele(stoppedVCs) )
////                mergedVC = ReferenceConfidenceVariantContextMerger.merge(stoppedVCs, gLoc, refBase, false, false, annotationEngine);
////            else
////                mergedVC = referenceBlockMerge(stoppedVCs, state, pos.getStart());
//			VariantContext variantContext = engine.variantCalling(values.iterator(),
//					parser.createGenomeLocation(chr, iter), genomeShare.getChromosomeInfo(chr));
//			if (variantContext == null)
//				continue;
//			CommonInfo info = variantContext.getCommonInfo();
//			HashMap<String, Object> maps = new HashMap<>();
//			maps.putAll(info.getAttributes());
//			maps.remove("SM");
//			info.setAttributes(maps);
//
//			outValue.set(variantContext, header);
//			context.write(NullWritable.get(), outValue);
//		}
	}
//	protected final class OverallState {
//        final LinkedList<VariantContext> VCs = new LinkedList<>();
//        final Set<String> samples = new HashSet<>();
//        GenomeLocation prevPos = null;
//        byte refAfterPrevPos;
//
//        public OverallState() {}
//    }
	
}
