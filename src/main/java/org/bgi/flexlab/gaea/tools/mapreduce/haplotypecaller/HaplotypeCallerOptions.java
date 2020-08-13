package org.bgi.flexlab.gaea.tools.mapreduce.haplotypecaller;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.PathFilter;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReferenceConfidenceMode;
import org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection.HaplotypeCallerArgumentCollection;
import org.bgi.flexlab.gaea.tools.haplotypecaller.pairhmm.PairHMM;
import org.seqdoop.hadoop_bam.SAMFormat;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class HaplotypeCallerOptions  extends GaeaOptions implements HadoopOptions{
	
	private final static String SOFTWARE_NAME = "HaplotypeCaller";
	
	private final static String SOFTWARE_VERSION = "1.0";
	
	private String region = null;
	
	private String reference = null;
	
	private int windowsSize = 10000;
	
	private int windowsExtends = 300;
	
	private int readShardSize = -1;
	
	private int readPaddingSize = 100;
	
	private List<String> userDisabledReadFilterNames = new ArrayList<>();
	
	private List<String> userEnabledReadFilterNames = new ArrayList<>();

	private  boolean disableToolDefaultReadFilters = false;

	private String dbsnp = null;
	
	private String allele = null;
	
	private boolean gvcfFormat = false;
	
	private HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
	
	private int reduceNumber = 100;
	
	private List<Path> inputs = new ArrayList<Path>();
	
	private String output = null;
	
	private boolean isSAM = false;
	
	private HashMap<String,String> comps  = new HashMap<String,String>();
	
	private int maxReadsPerPosition = 0;

	private boolean outputAllWindows;

	public List<Integer> GVCFGQBands = new ArrayList<>(70);
	
	public HaplotypeCallerOptions() {
		addOption("a","allSitePLs",false,"Annotate all sites with PLs");
		addOption("A","annotateNDA",false,"If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site");
		addOption("b","hets",true,"Heterozygosity value used to compute prior likelihoods for any locus");
		addOption("B","indel_hets",true,"Heterozygosity for indel calling");
		addOption("C","sample_ploidy",true,"Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).");
		addOption("c","shard_size",true,"read shard size.");
		addOption("d","shard_padding_size",true,"read shard padding size.");
		addOption("x","max_reads",true,"max reads for pileup.");
		addOption("E","windowExtendSize",true,"key window extend size.");
		addOption("e","max_depth_for_assembly",true,"max depth for assembly.");
		addOption("f", "format", false, "output format is gvcf");
		addOption("G", "gt_mode",true,"Specifies how to determine the alternate alleles to use for genotyping(DISCOVERY or GENOTYPE_GIVEN_ALLELES)");
		addOption("i", "input", true, "a bam or bam list for input", true);
		addOption("I","include_non_variant",false,"Include loci found to be non-variant after genotyping");
		addOption("j","heterozygosity_stdev",true,"Standard deviation of eterozygosity for SNP and indel calling");
		addOption("k", "knowSite", true, "known snp/indel file,the format is VCF4");
		addOption("K","dontIncreaseKmerSizes",false,"dont increase kmer sizes for cycles.");
		addOption("n", "reducer", true, "reducer numbers[100]");
		addOption("m","max_num_PL_values",true,"Maximum number of PL values to output");
		addOption("M","max_alternate_alleles",true,"Maximum number of alternate alleles to genotype");
		addOption("o", "output", true, "output directory", true);
		addOption("O","output_mode",true,"output mode(EMIT_VARIANTS_ONLY,EMIT_ALL_CONFIDENT_SITES,EMIT_ALL_SITES)");
		addOption("p","input_prior",true,"Input prior for calls(separation by Comma(,))");
		addOption("P","pairHMM",true,"pairHMM implementation:EXACT,ORIGINAL,LOGLESS_CACHING,AVX_LOGLESS_CACHING,AVX_LOGLESS_CACHING_OMP,FASTEST_AVAILABLE.");
		addOption("r", "reference", true, "reference index(generation by GaeaIndex) file path", true);
		addOption("R", "region", true, "One or more genomic intervals over which to operate");
		addOption("s","stand_emit_conf",true,"The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold");
		addOption("S","stand_call_conf",true,"The minimum phred-scaled confidence threshold at which variants should be called");
		addOption("u","uniquifySamples",false,"Assume duplicate samples are present and uniquify all names with '.variant' and file number index");
		addOption("U","useNewAFCalculator",false,"Use new AF model instead of the so-called exact model");
		addOption("w", "keyWindow", true, "window size for key[10000]");
		addOption("W", "outputAllWindow", false, "output N or uncover region windows [false]");
		FormatHelpInfo(SOFTWARE_NAME,SOFTWARE_VERSION);
	}
	
	private void initializeReadFilter() {
		
	}
	
	@Override
	public void setHadoopConf(String[] args, Configuration conf) {
		conf.setStrings("args", args);
	}

	@Override
	public void getOptionsFromHadoopConf(Configuration conf) {
		String[] args = conf.getStrings("args");
		this.parse(args);
	}

	@Override
	public void parse(String[] args) {
		try {
			cmdLine = parser.parse(options, args);
		} catch (ParseException e) {
			helpInfo.printHelp("Options:", options, true);
			System.exit(1);
		}

		if(args.length == 0 || getOptionBooleanValue("h", false)) {
			printHelpInfotmation(SOFTWARE_NAME);
			System.exit(1);
		}
		
		try {
			parseInput( getOptionValue("i",null));
		} catch (IOException e) {
			throw new UserException(e.toString());
		}

		if(inputs.size() <= 0)
			throw new RuntimeException("Input is empty!");

		if(getOptionBooleanValue("f",false)) {
			this.hcArgs.emitReferenceConfidence = ReferenceConfidenceMode.GVCF;
			this.gvcfFormat = true;
		}
		if(getOptionBooleanValue("K",false))
			this.hcArgs.assemblerArgs.dontIncreaseKmerSizesForCycles = true;
		this.hcArgs.maxDepthForAssembly = getOptionIntValue("e",0);
		
		this.windowsSize = getOptionIntValue("w",10000);
		this.windowsExtends = getOptionIntValue("E",300);
		this.reduceNumber = getOptionIntValue("n",100);
		this.readShardSize = getOptionIntValue("c",-1);
		this.readPaddingSize = getOptionIntValue("d",100);
		this.maxReadsPerPosition = getOptionIntValue("x",0);		
		this.output = getOptionValue("o",null);
		this.region = getOptionValue("R",null);
		this.reference = getOptionValue("r",null);
		this.dbsnp = getOptionValue("k",null);
		setOutputAllWindows(getOptionBooleanValue("W", false));
		
		if(dbsnp != null) {
			comps.put("DB", dbsnp);
		}
		
		setPairHMM(getOptionValue("P","AVX_LOGLESS_CACHING_OMP"));

	}
	
	private void parseInput(String input) throws IOException {
		Path path = new Path(input);
		FileSystem inFS = path.getFileSystem(new Configuration());
		PathFilter filter = file -> !file.getName().startsWith("_");
		if(inFS.isDirectory(path)){
			FileStatus[] stats = inFS.listStatus(path, filter);
			if(stats.length <= 0){
				System.err.println("Input File Path is empty! Please check input : " +path.toString());
				System.exit(-1);
			}
			for (FileStatus f: stats)
				inputs.add(f.getPath());
			return;
		}

		SAMFormat fmt = SAMFormat.inferFromData(inFS.open(path));
		
		if(fmt == SAMFormat.BAM)
			inputs.add(path);
		else {
			LineReader reader = new LineReader(inFS.open(path));
			Text line = new Text();
			
			while(reader.readLine(line) > 0 && line.getLength() != 0) {
				inputs.add(new Path(line.toString()));
			}
			reader.close();
		}
	}

	public String getRegion(){
		return this.region;
	}
	
	public String getReference(){
		return this.reference;
	}
	
	public int getWindowSize(){
		return windowsSize;
	}
	
	public int getWindowsExtendSize(){
		return windowsExtends;
	}
	
	public int getReadShardSize() {
		return this.readShardSize;
	}
	
	public int getReadShardPadding() {
		return this.readPaddingSize;
	}
	
	public List<String> getUserDisabledReadFilterNames(){
		return this.userDisabledReadFilterNames;
	}
	
	public List<String> getUserEnabledReadFilterNames(){
		return this.userEnabledReadFilterNames;
	}
	
	public boolean getDisableToolDefaultReadFilters() {
		return this.disableToolDefaultReadFilters;
	}
	
	public String getDBSnp() {
		return dbsnp;
	}
	
	public String getAlleleFile() {
		return this.allele;
	}
	
	public boolean isGVCF() {
		return this.gvcfFormat;
	}
	
	public int getReducerNumber() {
		return this.reduceNumber;
	}
	
	public List<Path> getInput(){
		return this.inputs;
	}
	
	public String getHeaderOutput() {
		return this.output;
	}
	
	public String getVCFOutput() {
		if(!output.endsWith("/"))
			output += "/";
		return output+"vcf";
	}
	
	public SAMFormat getInputFormat() {
		if(isSAM)
			return SAMFormat.SAM;
		return SAMFormat.BAM;
	}

	public boolean isOutputAllWindows() {
		return outputAllWindows;
	}

	public void setOutputAllWindows(boolean outputAllWindows) {
		this.outputAllWindows = outputAllWindows;
	}

	public HaplotypeCallerArgumentCollection getHaplotypeCallerArguments() {
		return this.hcArgs;
	}
	
	public List<String> getCompNames(){
		List<String> list = new ArrayList<String>();
		
		for(String name : comps.keySet())
			list.add(name);
		
		return list;
	}
	
	public int getMaxReadsPerPosition(){
		return this.maxReadsPerPosition;
	}
	
	private void setPairHMM(String args) {
		if(args == null)
			return;
		if(args.equals("EXACT"))
			hcArgs.likelihoodArgs.pairHMM = PairHMM.Implementation.EXACT;
		else if(args.equals("ORIGINAL"))
			hcArgs.likelihoodArgs.pairHMM = PairHMM.Implementation.ORIGINAL;
		else if(args.equals("ORIGINAL"))
			hcArgs.likelihoodArgs.pairHMM = PairHMM.Implementation.ORIGINAL;
		else if(args.equals("LOGLESS_CACHING"))
			hcArgs.likelihoodArgs.pairHMM = PairHMM.Implementation.LOGLESS_CACHING;
		else if(args.equals("AVX_LOGLESS_CACHING"))
			hcArgs.likelihoodArgs.pairHMM = PairHMM.Implementation.AVX_LOGLESS_CACHING;
		else if(args.equals("AVX_LOGLESS_CACHING_OMP"))
			hcArgs.likelihoodArgs.pairHMM = PairHMM.Implementation.AVX_LOGLESS_CACHING_OMP;
		else if(args.equals("FASTEST_AVAILABLE"))
			hcArgs.likelihoodArgs.pairHMM = PairHMM.Implementation.FASTEST_AVAILABLE;
		else
			throw new UserException.BadArgumentValueException("pairHMM",args);
	}
}
