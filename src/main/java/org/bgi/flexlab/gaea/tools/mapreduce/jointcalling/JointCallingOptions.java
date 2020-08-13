package org.bgi.flexlab.gaea.tools.mapreduce.jointcalling;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;
import org.bgi.flexlab.gaea.data.structure.vcf.AbstractVCFLoader;
import org.bgi.flexlab.gaea.tools.jointcalling.UnifiedGenotypingEngine.GenotypingOutputMode;
import org.bgi.flexlab.gaea.tools.jointcalling.UnifiedGenotypingEngine.OutputMode;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RealignerExtendOptions;
import org.seqdoop.hadoop_bam.VCFFormat;

import htsjdk.variant.vcf.VCFCodec;

public class JointCallingOptions extends GaeaOptions implements HadoopOptions{
	private final static String SOFTWARE_NAME = "JointCallingPrepare";
	private final static String SOFTWARE_VERSION = "1.0";
	
	private String inputList=null;
	
	private int samplePloidy = 2;//C
	private int MAX_ALTERNATE_ALLELES = 6;//M
	private int MAX_NUM_PL_VALUES = 100;//m
	private int windows_size = 5000;//w
	private int num_reducer = 100;//n
	private int num_mapper = 100;//N
	private boolean large_mode=false;//D

	private String tmpOut=null;//t
	private String output = null;//o
	private String reference = null;//r
	private String dbsnp = null;//k	
	private String MAGIC_HEADER_LINE = VCFCodec.VCF4_MAGIC_HEADER;//F
	
	//private String winFile=null;//W
	
	private double snpHeterozygosity = 1e-3;//b	
	private double indelHeterozygosity = 1.0/8000;//B
	
	private List<Path> input = new ArrayList<Path>();//i
	private List<Double> inputPrior = new ArrayList<Double>();//p
	private OutputMode outputMode = OutputMode.EMIT_VARIANTS_ONLY;//O
	private GenotypingOutputMode genotypeMode = GenotypingOutputMode.DISCOVERY;//G
	
	private boolean unquifySamples = false;//u
	private boolean bcfFormat = false;//f

	public boolean ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = false;//A
	public boolean ANNOTATE_ALL_SITES_WITH_PL = false;//a
	public boolean INCLUDE_NON_VARIANT = false;//I
	public boolean USE_NEW_AF_CALCULATOR = false; //U
	public double STANDARD_CONFIDENCE_FOR_CALLING = 30.0;//S
	public double STANDARD_CONFIDENCE_FOR_EMITTING = 30.0;//s
	public double heterozygosityStandardDeviation = 0.01;//j
	public String vcfHeaderFile = null;
	
	
	public JointCallingOptions(){
		addOption("a","allSitePLs",false,"Annotate all sites with PLs");
		addOption("A","annotateNDA",false,"If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site");
		addOption("b","hets",true,"Heterozygosity value used to compute prior likelihoods for any locus");
		addOption("B","indel_hets",true,"Heterozygosity for indel calling");
		addOption("C","sample_ploidy",true,"Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).");
		addOption("f", "format", false, "output format is bcf");
		//addOption("F", "vcfformat", true, "input vcf format version");
		addOption("G", "gt_mode",true,"Specifies how to determine the alternate alleles to use for genotyping(DISCOVERY or GENOTYPE_GIVEN_ALLELES)");
		//addOption("i", "input", true, "a gvcf or a gvcf list for input", true);
		addOption("i", "input", true, "a gvcf list for input", true);
		addOption("I","include_non_variant",false,"Include loci found to be non-variant after genotyping");
		addOption("j","heterozygosity_stdev",true,"Standard deviation of eterozygosity for SNP and indel calling");
		addOption("k", "knowSite", true, "known snp/indel file,the format is VCF4");
		addOption("D","largeSampleSize",false,"mode4 for large sample size(e.g. more than 1000 samples)");
		addOption("N", "mapper", true, "mapper number[100]");
		addOption("n", "reducer", true, "reducer numbers[100]");
		addOption("m","max_num_PL_values",true,"Maximum number of PL values to output");
		addOption("M","max_alternate_alleles",true,"Maximum number of alternate alleles to genotype");
		addOption("t","tmp_outdir",true,"directory for temporary files");
		addOption("o", "output", true, "output directory", true);
		addOption("O","output_mode",true,"output mode(EMIT_VARIANTS_ONLY,EMIT_ALL_CONFIDENT_SITES,EMIT_ALL_SITES)");
		addOption("p","input_prior",true,"Input prior for calls(separation by Comma(,))");
		addOption("r", "reference", true, "reference index(generation by GaeaIndex) file path", true);
		addOption("s","stand_emit_conf",true,"The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold");
		addOption("S","stand_call_conf",true,"The minimum phred-scaled confidence threshold at which variants should be called");
		addOption("u","uniquifySamples",false,"Assume duplicate samples are present and uniquify all names with '.variant' and file number index");
		addOption("U","useNewAFCalculator",false,"Use new AF model instead of the so-called exact model");
		addOption("w", "keyWindow", true, "window size for key[5000]");
		//addOption("W","windowsFile",true,"whole genome bed file"); //no need for window file input, it will be created by program
		FormatHelpInfo(SOFTWARE_NAME,SOFTWARE_VERSION);
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
			inputList=getOptionValue("i",null);
			//parseInput( getOptionValue("i",null));
			parseInput(inputList);
		} catch (IOException e) {
			throw new UserException(e.toString());
		}

		if(input.size() <= 0)
			throw new RuntimeException("Input is empty!");
		
		this.samplePloidy = getOptionIntValue("C",2);
		this.MAX_ALTERNATE_ALLELES = getOptionIntValue("M",6);
		this.MAX_NUM_PL_VALUES = getOptionIntValue("m",100);
		this.windows_size = getOptionIntValue("w",5000);
		this.num_reducer = getOptionIntValue("n",100);
		this.num_mapper=getOptionIntValue("N", 100);
		
		this.tmpOut=getOptionValue("t",null);
		this.output = getOptionValue("o",null);
		this.reference = getOptionValue("r",null);
		this.dbsnp = getOptionValue("k",null);
		parseInputPrior(getOptionValue("p",null));
		parseOutputMode(getOptionValue("O",OutputMode.EMIT_VARIANTS_ONLY.toString()));
		parseGenotypeMode(getOptionValue("G",GenotypingOutputMode.DISCOVERY.toString()));

		//this.winFile=getOptionValue("W",null);//added by gc. Windows file of ref
		this.large_mode=getOptionBooleanValue("D",false);
		this.ANNOTATE_ALL_SITES_WITH_PL = getOptionBooleanValue("a",false);
		this.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = getOptionBooleanValue("A",false);
		this.INCLUDE_NON_VARIANT = getOptionBooleanValue("I",false);
		this.unquifySamples = getOptionBooleanValue("u",false);
		this.USE_NEW_AF_CALCULATOR=getOptionBooleanValue("U",false);
		this.bcfFormat = getOptionBooleanValue("f",false);
		
		this.snpHeterozygosity = getOptionDoubleValue("b",1e-3);
		this.indelHeterozygosity = getOptionDoubleValue("B",1.0/8000);
		this.STANDARD_CONFIDENCE_FOR_EMITTING = getOptionDoubleValue("s",30.0);
		this.STANDARD_CONFIDENCE_FOR_CALLING = getOptionDoubleValue("S",30.0);
		this.heterozygosityStandardDeviation = getOptionDoubleValue("j",0.01);
	}
	
	private void parseOutputMode(String mode){
		if(mode.equals(OutputMode.EMIT_ALL_CONFIDENT_SITES.toString()))
			this.outputMode = OutputMode.EMIT_ALL_CONFIDENT_SITES;
		else if(mode.equals(OutputMode.EMIT_ALL_SITES.toString()))
			this.outputMode = OutputMode.EMIT_ALL_SITES;
		else if(mode.equals(OutputMode.EMIT_VARIANTS_ONLY.toString()))
			this.outputMode = OutputMode.EMIT_VARIANTS_ONLY;
		else
			throw new UserException.BadArgumentValueException("O",mode);
	}
	
	private void parseGenotypeMode(String mode){
		if(mode.equals(GenotypingOutputMode.DISCOVERY.toString()))
			this.genotypeMode = GenotypingOutputMode.DISCOVERY;
		else if(mode.equals(GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES.toString()))
			this.genotypeMode = GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES;
		else
			throw new UserException.BadArgumentValueException("G",mode);
	}
	
	private void parseInputPrior(String prior){
		if(prior == null || prior.length() == 0)
			return;
		String[] str = prior.trim().split(",");
		
		for(String s : str){
			inputPrior.add(Double.parseDouble(s));
		}
	}
	
	private void parseInput(String inputpath) throws IOException{
		if(inputpath == null)
			throw new UserException.BadArgumentValueException("i",inputpath);
		Path path = new Path(inputpath);
		Configuration conf = new Configuration();
		FileSystem inFS = path.getFileSystem(conf);
		PathFilter filter = file -> !file.getName().startsWith("_");
		if(inFS.isDirectory(path)){
			FileStatus[] fileStatuses = inFS.listStatus(path, filter);
			for (FileStatus f: fileStatuses) {
				if(f.getLen() <= 0)
					continue;
				input.add(f.getPath());
			}
			Path vcfHeaderPath = new Path(path.getParent().toString() + "/vcfFileHeader.vcf");
			if(inFS.exists(vcfHeaderPath))
				setVcfHeaderFile(vcfHeaderPath.toString());
			return;
		}

		boolean isvcf = AbstractVCFLoader.isVCFStream(inFS.open(path), MAGIC_HEADER_LINE);

		if(isvcf){
			input.add(path);
		}
		else{
			FSDataInputStream currInput;
			try {
				currInput = inFS.open(path);
				@SuppressWarnings("resource")
				LineReader lineReader = new LineReader(currInput,conf);
				Text line = new Text();
				
				while(lineReader.readLine(line) > 0){
					input.add(new Path(line.toString()));
				}
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}
			
			if(currInput != null)
				currInput.close();
		}
	}

	public int getSamplePloidy(){
		return this.samplePloidy;
	}
	
	public double getSNPHeterozygosity(){
		return this.snpHeterozygosity;
	}
	
	public double getINDELHeterozygosity(){
		return this.indelHeterozygosity;
	}
	
	public List<Double> getInputPrior(){
		return this.inputPrior;
	}
	
	public OutputMode getOutputMode(){
		return this.outputMode;
	}
	
	public GenotypingOutputMode getGenotypingOutputMode(){
		return this.genotypeMode;
	}
	
	public int getMaxAlternateAllele (){
		return this.MAX_ALTERNATE_ALLELES;
	}
	
	public int getMaxNumberPLValues(){
		return this.MAX_NUM_PL_VALUES;
	}
	
	public int getWindowsSize(){
		return this.windows_size;
	}
	public String getTmpOut() {
		return this.tmpOut;
	}
	public String getOutput(){
		if(output.endsWith("/"))
			return this.output+"vcf";
		else
			return this.output+"/vcf";
	}
	//added by gc
//	public String getWinFile() {
//		return this.winFile;
//	}
	public String getOutDir(){
		return this.output;
	}
	public List<Path> getInput(){
		return this.input;
	}
	
	public String getReference(){
		return this.reference;
	}
	
	public String getDBSnp(){
		return this.dbsnp;
	}
	
	public boolean isUniquifySamples(){
		return this.unquifySamples;
	}
	
	public int getReducerNumber(){
		return this.num_reducer;
	}
	public int getMapperNumber() {
		return this.num_mapper;
	}
	public boolean getMapperMode() {
		return this.large_mode;
	}
	public String getVcfHeaderFile() {
		return vcfHeaderFile;
	}

	public void setVcfHeaderFile(String vcfHeaderFile) {
		this.vcfHeaderFile = vcfHeaderFile;
	}

	public String getVCFHeaderOutput(){
		return this.output;
	}
	
	public String getInputList() {
		return this.inputList;
	}
	public void setInputList(String inputList){
		this.inputList=inputList;
	}
	public VCFFormat getOuptputFormat() {
	    if(bcfFormat)
	        return VCFFormat.BCF;
	    return VCFFormat.VCF;
	}
	public void setS(double v){
		this.STANDARD_CONFIDENCE_FOR_CALLING=v;
	}
	public void sets(double v){
		this.STANDARD_CONFIDENCE_FOR_EMITTING=v;
	}
}
