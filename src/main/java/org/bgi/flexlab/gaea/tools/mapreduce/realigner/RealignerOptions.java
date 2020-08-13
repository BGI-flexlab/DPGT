/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.mapreduce.realigner;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;
import org.seqdoop.hadoop_bam.SAMFormat;

public class RealignerOptions extends GaeaOptions implements HadoopOptions{
	private final static String SOFTWARE_NAME = "Realignment";
	private final static String SOFTWARE_VERSION = "1.0";
	
	private int winSize;
	private int reducerNumbers;
	private int maxReadsAtWindows;
	private int extendSize;
	private int snpWindowSize;
	private int minReadsAtPileup;
	private int maxIntervalSize;
	private int maxInsertSize;
	
	private String knowVariant;
	private String input;
	private String output;
	private String reference;
	
	private boolean samFormat;
	private boolean multiSample;
	
	private double mismatchThreshold = 0.0;
	private double LOD = 5.0;
	private AlternateConsensusModel model = AlternateConsensusModel.READS;
	
	public enum AlternateConsensusModel {
		/**
		 * generate alternate consensus model
		 */
		
		/**
		 * uses known indels only.
		 */
		DBSNP,
		/**
		 * uses know indels and the reads.
		 */
		READS,
		/**
		 * uses know indels and the reads and 'Smith-Waterman'.
		 */
		SW
	}
	
	public RealignerOptions(){
		addOption("c", "consensusModel",true,"Determines how to compute the possible alternate consenses.model:DBSNP,READS.[READS]");
		addOption("d", "LOD",true,"LOD threshold above which the cleaner will clean [5.0].");
		addOption("e", "windowExtendSize", true, "window extend size[500]");
		addOption("i", "input", true, "input directory", true);
		addOption("I", "insertSize",true,"maximum insert size of read pairs that we attempt to realign [3000].");
		addOption("k", "knowSite", true, "known snp/indel file,the format is VCF4");
		addOption("l", "minReads", true, "minimum reads at a locus to enable using the entropy calculation[4].");
		addOption("L", "intervalLength", true, "max interval length[500].");
		addOption("M", "multiSample", false, "mutiple sample realignment[false]");
		addOption("m", "maxReadsAtWindows", true, "max reads numbers at on windows[1000000]");
		addOption("n", "reducer", true, "reducer numbers[30]");
		addOption("o", "output", true, "output directory", true);
		addOption("r", "reference", true, "reference index(generation by GaeaIndex) file path", true);
		addOption("s", "samformat", false, "input file is sam format");		
		addOption("t", "mismatch", true, "fraction of base qualities needing to mismatch for a position to have high entropy[0]");
		addOption("w", "keyWindow", true, "window size for key[10000]");
		addOption("W", "window", true, "window size for calculating entropy or SNP clusters[10]");
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
			System.err.println(e.getMessage());
			printHelpInfotmation(RealignerExtendOptions.SOFTWARE_NAME);
			System.exit(1);
		}
		
		input = getOptionValue("i",null);
		output = getOptionValue("o",null);
		reference = getOptionValue("r",null);
		knowVariant = getOptionValue("k",null);
		
		if(!getOptionValue("c","READS").equals("READS"))
			model = AlternateConsensusModel.DBSNP;
		
		winSize = getOptionIntValue("w",10000);
		reducerNumbers = getOptionIntValue("n",30);
		maxReadsAtWindows = getOptionIntValue("m",1000000);
		extendSize = getOptionIntValue("e",500);
		snpWindowSize = getOptionIntValue("W",10);
		minReadsAtPileup = getOptionIntValue("l",4);
		maxIntervalSize = getOptionIntValue("L",500);
		maxInsertSize = getOptionIntValue("I",3000);
		
		mismatchThreshold = getOptionDoubleValue("t",0);
		LOD = getOptionDoubleValue("d",5.0);
		
		samFormat = getOptionBooleanValue("s",false);
		multiSample = getOptionBooleanValue("M",false);
		
		if(output != null && !output.endsWith("/")){
			output += "/";
		}
	}

	public boolean isMultiSample(){
		return multiSample;
	}
	
	public SAMFormat getInputFormat(){
		if(samFormat)
			return SAMFormat.SAM;
		return SAMFormat.BAM;
	}
	
	public int getWindowsSize(){
		return winSize;
	}
	
	public int getReducerNumber(){
		return reducerNumbers;
	}
	
	public int getMaxReadsAtWindows(){
		return this.maxReadsAtWindows;
	}
	
	public String getReference(){
		return reference;
	}
	
	public String getKnowVariant(){
		return knowVariant;
	}
	
	public String getRealignerInput(){
		return input;
	}
	
	public String getRealignerOutput(){
		return output+"realigner";
	}
	
	public String getFixmateInput(boolean control){
		if(control)
			return getRealignerOutput();
		else
			return getRealignerInput();
	}
	
	public String getFixmateOutput(){
		return output+"fixmate";
	}
	
	public int getExtendSize(){
		return this.extendSize;
	}
	
	public int getSNPWindowSize(){
		return this.snpWindowSize;
	}
	
	public double getMismatchThreshold(){
		return mismatchThreshold;
	}
	
	public int getMinReads(){
		return minReadsAtPileup;
	}
	
	public int getMaxInterval(){
		return maxIntervalSize;
	}
	
	public double getLODThreshold(){
		return this.LOD;
	}
	
	public AlternateConsensusModel getConsensusModel(){
		return model;
	}
	
	public int getMaxInsertSize(){
		return this.maxInsertSize;
	}
}
