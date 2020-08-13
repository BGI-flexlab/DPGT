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
package org.bgi.flexlab.gaea.tools.mapreduce.fastqqualitycontrol;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;

public class FastqQualityControlOptions extends GaeaOptions implements
		HadoopOptions {
	private final static String SOFTWARE_NAME = "FastqQualityControl";
	private final static String SOFTWARE_VERSION = "1.0";
	private final static int LENGTH_QUALITY = 64;

	public FastqQualityControlOptions() {
		addOption("1", "in1fq", true, "the first fastq");
		addOption("2", "in2fq", true, "the second fastq");
		addOption("3", "adapter1", true, "the first adapter list file");
		addOption("4", "adapter2", true, "the second adapter list file");
		addOption("l", "lowQuality", true, "low quality value (default:5)");
		addOption(
				"Q",
				"qualitySystem",
				true,
				"Quality score schema(default:1)\n0: NCBI/Sanger or Illumina 1.8 and later; "
						+ " \n1: Illumina Pipeline 1.5 to 1.7;  \n2: Illumina Pipeline 1.3 and 1.4;  "
						+ "\n3: Illumina Pipeline 1.2 and earlier");
		addOption("C", "sanger", false, "change base quality to sanger");
		addOption("q", "qualityRate", true, "low quality rate(default:0.5)");
		addOption("t", "readType", true, "read name type.(default:0)\n0: reads_xxx/1\n1: reads_xx: 1:N:XX\n2: reads_xx");
		addOption("N", "NRate", true, "Maximum N rate(default:0.1)");
		addOption("s", "trim", true,
				"trim some bp of the read's head and tail, they means: (read1's head and tail and read2's head and tail  [0,0,0,0]");
		addOption("5", "ignored1fq", false, "not output reads from fastq1");
		addOption("6", "ignored2fq", false, "not output reads from fastq2");
		addOption("m", "multiSample", true, "Mulit samples list");
		addOption("n", "reducerNumber", true, "reducer number");
		addOption("o", "output", true, "output directory", true);
		addOption("d", "qualityFreq", false, "output quality frequency statice");
		addOption("M", "ignoredSample", false, "QC statictis together");
		addOption("y", "dyncut", false, "run dyncutadaptor process");
		addOption("S", "seed", true, "initial length of adaptor.(default:10)");
		addOption("T", "ignoredTail", true,
				"don't cut tail in the last n bp.(default:3)");
		addOption("A", "mismatch", true, "tolerate mismatchs.(default:1)");
		addOption("B", "minimum", true, "don't keep seqences shorter.(default:0)");
		addOption("h", "help", false, "help information");

		FormatHelpInfo(SOFTWARE_NAME, SOFTWARE_VERSION);
	}

	private String adapter1;
	private String adapter2;
	private String input1Fastq;
	private String input2Fastq;
	private String outputDir;
	private String multiSampleList;

	private double qualRate;
	private double NRate;

	private int minimumQualityScore;
	private int Q;
	private int[] trim;
	private int reducerNumber;
	private int lowQual;
	private int cest;
	private int minimum;
	private int bias;
	private int length_adaptor;
	private int readType;

	private boolean SEdata = false;
	private boolean filterSE;
	private boolean qualTrim;
	private boolean dyncut = false;
	private boolean multiStatis = true;
  private boolean qualFreq = false;
	private boolean ignoredfastq1;
	private boolean ignoredfastq2;

	@Override
	public void parse(String[] args) {
		try {
			cmdLine = parser.parse(options, args);
		} catch (ParseException e) {
			System.err.println(e.getMessage());
			printHelpInfotmation(SOFTWARE_NAME);
			System.exit(1);
		}
		
		if(!options.hasOption("1") && !options.hasOption("m")){
			System.err.println("must set -1 or -m");
			System.exit(1);
		}

		input1Fastq = getOptionValue("1", null);
		input2Fastq = getOptionValue("2", null);
		adapter1 = getOptionValue("3", null);
		adapter2 = getOptionValue("4", null);
		outputDir = getOptionValue("o", null);
		multiSampleList = getOptionValue("m", null);

		qualRate = getOptionDoubleValue("q", 0.5);
		NRate = getOptionDoubleValue("N", 0.1);

		Q = getOptionIntValue("Q", 1);
		setTrim(getOptionValue("s", "0,0,0,0"));
		reducerNumber = getOptionIntValue("n", 100);
		lowQual = getOptionIntValue("l", 5);
		cest = getOptionIntValue("T", 3);
		minimum = getOptionIntValue("B", 30);
		bias = getOptionIntValue("A", 1);
		length_adaptor = getOptionIntValue("S", 10);
		readType = getOptionIntValue("t",0);

		qualTrim = getOptionBooleanValue("C", false);
		dyncut = getOptionBooleanValue("y", false);
		multiStatis = getOptionBooleanValue("M", true);
		qualFreq = getOptionBooleanValue("d", false);
		ignoredfastq1 = getOptionBooleanValue("5", false);
		ignoredfastq2 = getOptionBooleanValue("6", false);

		setQualitySystem(Q);
	}

	private void setQualitySystem(int Q) {
		if (Q == 0) { // NCBI/Sanger or Illumina 1.8 and later
			minimumQualityScore = 33;
			System.out
					.println("> Quality score schema: NCBI/Sanger or Illumina 1.8 and later (ASCII 33 to 93)");
		} else if (Q == 1) { // Illumina Pipeline 1.5 to 1.7
			minimumQualityScore = 64;
			System.out
					.println("> Quality score schema: Illumina Pipeline 1.5 to 1.7 (ASCII 64 to 104)");
		} else if (Q == 2) { // Illumina Pipeline 1.3 and 1.4
			minimumQualityScore = 64;
			System.out
					.println("> Quality score schema: Illumina Pipeline 1.3 and 1.4 (ASCII 64 to 104, Values 0 (@) and 1 (A) are not used anymore. Value 2 (B) has special meaning and is used as a trim clipping.)");
		} else if (Q == 3) { // Illumina Pipeline 1.2 and earlier
			minimumQualityScore = 59;
			System.out
					.println("> Quality score schema: Illumina Pipeline 1.2 and earlier (ASCII 59 to 104)");
		} else {
			System.out
					.println("> Invailid value of quality score schema. Please check it. (0: NCBI/Sanger or Illumina 1.8 and later; 1: Illumina Pipeline 1.5 to 1.7; 2: Illumina Pipeline 1.3 and 1.4; 3: Illumina Pipeline 1.2 and earlier.)");
			System.exit(1);
		}
	}

	public String getAdapter1() {
		return adapter1;
	}

	public String getAdapter2() {
		return adapter2;
	}

	public String getInputFastq1() {
		return input1Fastq;
	}

	public String getInputFastq2() {
		return input2Fastq;
	}

	public String getOutputDirectory() {
		return outputDir;
	}

	public int getLowQuality() {
		return lowQual;
	}

	public double getQualityRate() {
		return qualRate;
	}

	public double getNRate() {
		return NRate;
	}

	public int getReducerNumber() {
		return reducerNumber;
	}

	public int getQuality() {
		return minimumQualityScore;
	}

	public int getQualityLength() {
		return LENGTH_QUALITY;
	}

	public void setQualityScore(int qual) {
		this.minimumQualityScore = qual;
	}

	public boolean isIgnoredFastq1() {
		return ignoredfastq1;
	}

	public boolean isIgnoredFastq2() {
		return ignoredfastq2;
	}

	public String getMultiSampleList() {
		return multiSampleList;
	}

	public boolean isSEdata() {
		return SEdata;
	}

	public void setSEdata(boolean sEdata) {
		SEdata = sEdata;
	}

	// trimStr : read1HeadTrim,read1TailTrim,read2HeadTrim,read2TailTrim  (int,int,int,int)
	private void setTrim(String trimStr) {
		trim = new int[4];
		String[] trimTemp = trimStr.trim().split(",");
		for (int i = 0; i < 4 ; i++){
			if(i < trimTemp.length) {
				trim[i] = Integer.parseInt(trimTemp[i]);
			}else {
				trim[i] = 0;
			}
		}
	}

	public int[] getTrim() {
		return trim;
	}

	public boolean isFilterSE() {
		return filterSE;
	}

	public boolean isDyncut() {
		return dyncut;
	}

	public int getCest() {
		return cest;
	}

	public int getBias() {
		return bias;
	}

	public int getAdaptorLength() {
		return length_adaptor;
	}

	public int getMinimum() {
		return minimum;
	}

	public boolean isQualityTrim() {
		return qualTrim;
	}
	
	public boolean isMultiStatis() {
		return this.multiStatis;
	}

	public int getQualitySystem() {
		return this.Q;
	}

  public boolean isQualityFrequency() {
		return qualFreq;
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
	
	public int getReadType(){
		return readType;
	}
}