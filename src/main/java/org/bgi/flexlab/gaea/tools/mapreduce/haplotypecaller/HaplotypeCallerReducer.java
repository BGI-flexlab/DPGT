package org.bgi.flexlab.gaea.tools.mapreduce.haplotypecaller;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.input.bed.RegionHdfsParser;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVariantContextWriter;
import org.bgi.flexlab.gaea.data.mapreduce.util.VariantContextHadoopWriter;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.data.structure.dbsnp.DbsnpShare;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.variant.filter.VariantRegionFilter;
import org.bgi.flexlab.gaea.tools.haplotypecaller.HaplotypeCallerTraversal;
import org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection.HaplotypeCallerArgumentCollection;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.haplotypecaller.writer.GVCFHadoopWriter;
import org.bgi.flexlab.gaea.util.Window;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;

public class HaplotypeCallerReducer extends Reducer<WindowsBasedWritable, SamRecordWritable, NullWritable, VariantContextWritable>{

	/**
     * region
     */
    private RegionHdfsParser region = null;
    
    private HaplotypeCallerOptions options = new HaplotypeCallerOptions();
    
    /**
     * sam file header
     */
    private SAMFileHeader header;

    /**
     * shared reference
     */
    private ReferenceShare genomeShare;
    
    /**
     * shared dbsnp
     */
    private DbsnpShare dbsnpShare = null;
    
    /**
     * vcf loader
     */
	private VCFLocalLoader DBloader = null;
	
	private DbsnpShare alleleShare = null;
	
	private VCFLocalLoader alleleLoader = null;
	
	/**
	 * vcf filter
	 */
	private VariantRegionFilter filter = null;
	
	private HaplotypeCallerTraversal haplotypecaller = null;
	private Map<String,HaplotypeCallerTraversal> sampleHaplotypecallers = null;
	private Map<String,GVCFHadoopWriter> gvcfWriters = null;

	/**
	 * variant context writer
	 */
	private GaeaVariantContextWriter writer = null;

	protected HashMap<Integer, String> sampleIDs = null;
    
	@Override
    protected void setup(Context context) throws IOException {
		Configuration conf = context.getConfiguration();
		options.getOptionsFromHadoopConf(conf);
		sampleHaplotypecallers = new HashMap<>();
		gvcfWriters = new HashMap<>();

		if(options.getRegion() != null){
			region = new RegionHdfsParser();
            region.parseBedFileFromHDFS(options.getRegion(), false);
		}
		
		header = SamHdfsFileHeader.getHeader(conf);
		List<SAMReadGroupRecord> list = header.getReadGroups();
		sampleIDs = new HashMap<>();
		for (int i = 0; i < list.size(); i++) {
			sampleIDs.put(i, list.get(i).getSample());
		}
        genomeShare = new ReferenceShare();
        genomeShare.loadChromosomeList(options.getReference());

        filter = new VariantRegionFilter();
        if(options.getDBSnp() != null) {
        	dbsnpShare = new DbsnpShare(options.getDBSnp(), options.getReference());
        	dbsnpShare.loadChromosomeList(options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
        	DBloader = new VCFLocalLoader(options.getDBSnp());
        }
        
        if(options.getAlleleFile() != null) {
        	alleleShare = new DbsnpShare(options.getAlleleFile(), options.getReference());
        	alleleShare.loadChromosomeList(options.getAlleleFile() + VcfIndex.INDEX_SUFFIX);
        	alleleLoader = new VCFLocalLoader(options.getAlleleFile());
        }

		HaplotypeCallerArgumentCollection hcArgs = options.getHaplotypeCallerArguments();
		for(String sample: sampleIDs.values()) {
        	if(!sampleHaplotypecallers.containsKey(sample)) {
//				hcArgs.sampleNameToUse = sample;
				SAMFileHeader sampleHeader = SamHdfsFileHeader.createHeaderFromSampleName(header, sample);
				haplotypecaller = new HaplotypeCallerTraversal(region, options, sampleHeader);
				sampleHaplotypecallers.put(sample, haplotypecaller);
			}
		}

		if (options.isGVCF()) {
			try {
				for(String sample: sampleIDs.values()) {
					if(!gvcfWriters.containsKey(sample)) {
						GVCFHadoopWriter writer = new GVCFHadoopWriter(context, haplotypecaller.getVCFHeader("multiSample"), hcArgs.GVCFGQBands, hcArgs.samplePloidy, sample);
						gvcfWriters.put(sample, writer);
					}
				}
//				writer = new GVCFHadoopWriter(context, haplotypecaller.getVCFHeader(), hcArgs.GVCFGQBands, hcArgs.samplePloidy);
			} catch (IllegalArgumentException e) {
				throw new UserException.BadArgumentValueException("GQBands", "are malformed: " + e.getMessage());
			}
		}else
			writer = new VariantContextHadoopWriter(context,haplotypecaller.getVCFHeader());

	}
	
	private ArrayList<VariantContext> getRegionVatiantContext(String chr,int number,int winSize,int end,DbsnpShare dbsnpShare,VCFLocalLoader loader){
		ArrayList<VariantContext> dbsnps = null;
		if(dbsnpShare != null) {
			long startPosition = dbsnpShare.getStartPosition(chr, number, winSize);
			if(startPosition >= 0)
				dbsnps = filter.loadFilter(loader, chr, startPosition, end);
		}
		
		return dbsnps;
	}
	
	private RefMetaDataTracker createTracker(String chr,int number,int winSize,int end) {
		RefMetaDataTracker tracker = null;
		
		if(dbsnpShare != null) {
			if(tracker == null)
				tracker = new RefMetaDataTracker();
			ArrayList<VariantContext> dbsnps = getRegionVatiantContext(chr,number,winSize,end,dbsnpShare,DBloader);
			if(dbsnps == null)
				dbsnps = new ArrayList<VariantContext>();
			tracker.add(RefMetaDataTracker.DB_VALUE, dbsnps);
		}
		if(alleleShare != null) {
			if(tracker == null)
				tracker = new RefMetaDataTracker();
			ArrayList<VariantContext> dbsnps = getRegionVatiantContext(chr,number,winSize,end,alleleShare,alleleLoader);
			if(dbsnps == null)
				dbsnps = new ArrayList<VariantContext>();
			tracker.add(RefMetaDataTracker.ALLELE_VALUE, dbsnps);
		}
		
		return tracker;
	}
	
	@Override
    public void reduce(WindowsBasedWritable key, Iterable<SamRecordWritable> values, Context context) throws IOException, InterruptedException {
		int index = key.getChromosomeIndex();
		if(index < 0)
			return;

		haplotypecaller = sampleHaplotypecallers.get(sampleIDs.get(key.getSampleID()));

		int start = key.getWindowsNumber() * options.getWindowSize();
		int contigLength = header.getSequenceDictionary().getSequence(index).getSequenceLength();
		int end = Math.min(contigLength, start + options.getWindowSize());
		Window win = new Window(header.getSequence(index).getSequenceName(), key.getChromosomeIndex(), start, end);
		String chr = win.getContigName();
		ChromosomeInformationShare chrInfo = genomeShare.getChromosomeInfo(chr,true);
		
		RefMetaDataTracker tracker = createTracker(chr,key.getWindowsNumber(),options.getWindowSize(),end);
		haplotypecaller.dataSourceReset(win, values, chrInfo, tracker);
		if(options.isGVCF()) {
			writer = gvcfWriters.get(sampleIDs.get(key.getSampleID()));
			haplotypecaller.traverse(writer);
			((GVCFHadoopWriter) writer).emitCurrentBlock();
		}else {
			haplotypecaller.traverse(writer);
		}
	}
	
	@Override
    protected void cleanup(Context context) {
		if(options.isGVCF()) {
			for(String sample: sampleIDs.values()) {
				if(gvcfWriters.containsKey(sample)) {
					gvcfWriters.get(sample).close();
				}
			}
		}else
			writer.close();
		haplotypecaller.clear();
    }
}
