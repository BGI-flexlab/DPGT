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
package org.bgi.flexlab.gaea.tools.mapreduce.bamqualitycontrol;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.PathFilter;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.GenericOptionsParser;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;
import org.seqdoop.hadoop_bam.SAMFormat;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class BamQualityControlOptions extends GaeaOptions implements HadoopOptions {
	private final static String SOFTWARE_NAME = "Bam Quality Control";
	private final static String SOFTWARE_VERSION = "1.0";
	
	public BamQualityControlOptions() {
		// TODO Auto-generated constructor stub
		addOption("i", "input", true, "Path of the alignment result file(s) or bam list", true);
		addOption("o", "output", true, "Path of the bamqc report file", true);
		addOption("d", "index", true, "Path of the reference index list", true);
		addOption("f", "ref", true, "reference file for read cram");
		addBooleanOption("b", "basic", "Basic mapping Statitis of each sample without coverage info");
		addBooleanOption("l", "lane", "Get basic information by lane rather than samples");
		addOption("r", "reducer", true, "The reducer nums");
		addOption("R", "region", true, "Region string(chr:start-end)");
		addOption("B", "bedFile", true, "bed format file");
		addBooleanOption("D", "distributeCache", "distribute cache the reference");
		addOption("a", "anRegion", true, "annotation regions");
//		addOption("s", "statRegion", true, "output detail infomation for this bedfile(WGS mode)");
		addOption("c", "cnvRegion", true, "cnv regions in bed format");
		addOption("n", "minDepth", true, "minimum depth for single region statistics");
		addOption("I", "insertSize", true, "insert size");
		addOption("P", "insertSizeWithoutDup", true, "insert size without duplication");
		addBooleanOption("C", "cnvDepth", "output the cnv depth file");
		addBooleanOption("x", "gender", "inhibit gender prediction");
		addBooleanOption("m", "unmapped", "inhibit unmapped region output");
		addOption("h", "help", false, "help information");
	}

	private String alignmentFilePath;

	private List<Path> inputs = new ArrayList<>();

	/**
	 * 一致性序列文件名称
	 */
	private String outputPath;

	/**
	 * 参考基因组序列文件路径
	 */
	private String referenceSequencePath;

	private String localReferenceSequencePath;

	/**
	 * Reducer个数
	 */
	private int reducerNum;
	
	/**
	 * region in string
	 */
	private String region;
	
	/**
	 * regions in bed file
	 */
	private String bedfile;
	
	/**
	 * basic?
	 */
	private boolean isBasic;
	
	/**
	 * lane?sample?
	 */
	private boolean islane;
	
	/**
	 * distribute cache
	 */
	private boolean distributeCache;
	
	/**
	 * single region
	 */
	private String singleRegion;

	/**
	 * minimum depth for singleRegion
	 */
	private int minSingleRegionDepth;
	
	/**
	 * do cnv depth ? 
	 */
	private boolean cnvDepth;
	
	/**
	 * gender predict;
	 */
	private boolean genderPredict;
	
	private boolean outputUnmapped;
	
	private String cnvRegion;
	
	private int intsertSize;
	
	private int insertSizeWithoutDup;

	@Override
	public void setHadoopConf(String[] args, Configuration conf) {
		// TODO Auto-generated method stub
		String[] otherArgs;
		try {
			otherArgs = new GenericOptionsParser(args).getRemainingArgs();
			conf.setStrings("args", otherArgs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@Override
	public void getOptionsFromHadoopConf(Configuration conf) {
		// TODO Auto-generated method stub
		String[] args = conf.getStrings("args");
		this.parse(args);
	}

	@Override
	public void parse(String[] args) {
		// TODO Auto-generated method stub
		try {
			cmdLine = parser.parse(options, args);
		} catch (ParseException e) {
			System.err.println(e.getMessage());
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

		cnvDepth = getOptionBooleanValue("C", false);
		cnvRegion = getOptionValue("c", null);
		outputUnmapped = getOptionBooleanValue("m", true);
		genderPredict = getOptionBooleanValue("x", true);
		singleRegion = getOptionValue("a", null);
		minSingleRegionDepth = getOptionIntValue("n", 0);
		distributeCache = getOptionBooleanValue("D", false);
		isBasic = getOptionBooleanValue("b", false);
		islane = getOptionBooleanValue("l", false);
		referenceSequencePath = getOptionValue("d", null);
		localReferenceSequencePath = getOptionValue("f", null);
		outputPath = getOptionValue("o", null);
		reducerNum = getOptionIntValue("r", 30);
		region = getOptionValue("R", null);
		bedfile = getOptionValue("B", null);
		intsertSize = getOptionIntValue("I", 2000);
		insertSizeWithoutDup = getOptionIntValue("P", 2000);
		if((cnvRegion == null) && cnvDepth == true) {
			System.err.println("cnv Region not setted, use -B to replace.");
			if(bedfile == "") {
				System.err.println("region don't exists.");
				System.exit(1);
			} else {
				cnvRegion = bedfile;
			}
		}
	}

	private void parseInput(String input) throws IOException {
		Path path = new Path(input);
		FileSystem inFS = path.getFileSystem(new Configuration());
		if(inFS.isDirectory(path)){
			PathFilter filter = file -> !file.getName().startsWith("_");
			FileStatus[] stats=inFS.listStatus(path, filter);
			if(stats.length <= 0){
				System.err.println("Input File Path is empty! Please check input : " +path.toString());
				System.exit(-1);
			}

			for (FileStatus f: stats)
				inputs.add(f.getPath());
			return;
		}

		SAMFormat fmt = SAMFormat.inferFromData(inFS.open(path));

		if(fmt == SAMFormat.BAM || input.endsWith(".cram")) {
			inputs.add(path);
		}
		else {
			LineReader reader = new LineReader(inFS.open(path));
			Text line = new Text();

			while(reader.readLine(line) > 0 && line.getLength() != 0) {
				inputs.add(new Path(line.toString()));
			}
			reader.close();
		}
	}

	public List<Path> getInputs() {
		return inputs;
	}

	public void setInputs(List<Path> inputs) {
		this.inputs = inputs;
	}

	/**
	 * @return the outputPath
	 */
	public String getOutputPath() {
		return outputPath;
	}

	public String getTempPath() {
		return this.outputPath + "/tmp";
	}
	
	/**
	 * @param outputPath the outputPath to set
	 */
	public void setOutputPath(String outputPath) {
		this.outputPath = outputPath;
	}

	/**
	 * @return the referenceSequencePath
	 */
	public String getReferenceSequencePath() {
		return referenceSequencePath;
	}

	/**
	 * @param referenceSequencePath the referenceSequencePath to set
	 */
	public void setReferenceSequencePath(String referenceSequencePath) {
		this.referenceSequencePath = referenceSequencePath;
	}

	public String getLocalReferenceSequencePath() {
		return localReferenceSequencePath;
	}

	/**
	 * @return the region
	 */
	public String getRegion() {
		return region;
	}

	/**
	 * @param region the region to set
	 */
	public void setRegion(String region) {
		this.region = region;
	}

	/**
	 * @return the bedfile
	 */
	public String getBedfile() {
		return bedfile;
	}

	/**
	 * @param bedfile the bedfile to set
	 */
	public void setBedfile(String bedfile) {
		this.bedfile = bedfile;
	}

	public int getReducerNum() {
		return reducerNum;
	}

	public void setReducerNum(int reducerNum) {
		this.reducerNum = reducerNum;
	}

	/**
	 * @return the isBasic
	 */
	public boolean isBasic() {
		return isBasic;
	}

	/**
	 * @param isBasic the isBasic to set
	 */
	public void setBasic(boolean isBasic) {
		this.isBasic = isBasic;
	}

	/**
	 * @return the islane
	 */
	public boolean isIslane() {
		return islane;
	}

	/**
	 * @param islane the islane to set
	 */
	public void setIslane(boolean islane) {
		this.islane = islane;
	}

	/**
	 * @return the distributeCache
	 */
	public boolean isDistributeCache() {
		return distributeCache;
	}

	/**
	 * @param distributeCache the distributeCache to set
	 */
	public void setDistributeCache(boolean distributeCache) {
		this.distributeCache = distributeCache;
	}

	/**
	 * @return the singleRegion
	 */
	public String getSingleRegion() {
		return singleRegion;
	}

	/**
	 * @param singleRegion the singleRegion to set
	 */
	public void setSingleRegion(String singleRegion) {
		this.singleRegion = singleRegion;
	}

	/**
	 * @return the minSingleRegionDepth
	 */
	public int getMinSingleRegionDepth() {
		return minSingleRegionDepth;
	}

	/**
	 * @param minSingleRegionDepth the minSingleRegionDepth to set
	 */
	public void setMinSingleRegionDepth(int minSingleRegionDepth) {
		this.minSingleRegionDepth = minSingleRegionDepth;
	}

	/**
	 * @return the cnvDepth
	 */
	public boolean isCnvDepth() {
		return cnvDepth;
	}

	/**
	 * @param cnvDepth the cnvDepth to set
	 */
	public void setCnvDepth(boolean cnvDepth) {
		this.cnvDepth = cnvDepth;
	}

	/**
	 * @return the genderPredict
	 */
	public boolean isGenderPredict() {
		return genderPredict;
	}

	/**
	 * @param genderPredict the genderPredict to set
	 */
	public void setGenderPredict(boolean genderPredict) {
		this.genderPredict = genderPredict;
	}

	/**
	 * @return the outputUnmapped
	 */
	public boolean isOutputUnmapped() {
		return outputUnmapped;
	}

	/**
	 * @param outputUnmapped the outputUnmapped to set
	 */
	public void setOutputUnmapped(boolean outputUnmapped) {
		this.outputUnmapped = outputUnmapped;
	}

	/**
	 * @return the cnvRegion
	 */
	public String getCnvRegion() {
		return cnvRegion;
	}

	/**
	 * @param cnvRegion the cnvRegion to set
	 */
	public void setCnvRegion(String cnvRegion) {
		this.cnvRegion = cnvRegion;
	}
	
	public int getInsertSzie() {
		return intsertSize;
	}
	
	public int getInsertSzieWithoutDup() {
		return insertSizeWithoutDup;
	}
	
	public boolean isGenderDepth() {
		return getSingleRegion() != null || 
				(getBedfile() != null || getRegion() != null);
	}
}
