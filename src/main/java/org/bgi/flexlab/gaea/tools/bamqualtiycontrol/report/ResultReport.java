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
package org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report;

import htsjdk.samtools.SAMRecordIterator;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.CompressionCodecFactory;
import org.apache.hadoop.io.compress.CompressionOutputStream;
import org.apache.hadoop.mapreduce.Reducer.Context;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.structure.positioninformation.depth.PositionDepth;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegion;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegion.Regiondata;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegionStatistic;
import org.bgi.flexlab.gaea.tools.mapreduce.bamqualitycontrol.BamQualityControlOptions;
import org.bgi.flexlab.gaea.util.SamRecordDatum;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;

public abstract class ResultReport {
	
	private SingleRegion singleRegion;
	
	protected BasicReport basicReport;
	
	protected BamQualityControlOptions options;
	
	protected RegionCoverReport regionCoverReport;
	
	protected RegionCoverReport rmdupRegionCoverReport;
	
	protected SingleRegionReport cnvSingleRegionReport;
	
	protected UnmappedReport unmappedReport;
	
	protected int[] insertSize;
	
	protected int[] insertSizeWithoutDup;
	
	protected ReferenceShare genome;
	
	public ResultReport(BamQualityControlOptions options) throws IOException {
		this.options = options;
		if(options.getSingleRegion() != null) {
			singleRegion = new SingleRegion();
			singleRegion.parseRegionsFileFromHDFS(options.getSingleRegion(), false, 0);
		}
		loadReference();
	}
	
	private void loadReference() {
		genome = new ReferenceShare();
		if(options.isDistributeCache()) {
			genome.loadChromosomeList();
		} else {
			genome.loadChromosomeList(options.getReferenceSequencePath());
		}
	}
	
	public void initReports() throws IOException {
		basicReport = new BasicReport();
		unmappedReport = new UnmappedReport();
		regionCoverReport = new RegionCoverReport(1000);
		rmdupRegionCoverReport = new RegionCoverReport(1000);
		
		if(options.getSingleRegion() != null) {
			cnvSingleRegionReport = new SingleRegionReport(singleRegion);
		}
		insertSize = new int[options.getInsertSzie()];
		insertSizeWithoutDup = new int[options.getInsertSzieWithoutDup()];
		Arrays.fill(insertSize, 0);
		Arrays.fill(insertSizeWithoutDup, 0);
	}
	
	public boolean unmappedReport(long winNum, String chrName, Iterable<Text> values) {
		return unmappedReport.constructMapReport(winNum, chrName, values, basicReport);
	}

	public boolean unmappedReport(long winNum, String chrName, SAMRecordIterator values) {
		return unmappedReport.constructMapReport(winNum, chrName, values, basicReport);
	}
	
	public boolean mappedReport(SamRecordDatum datum, String chrName, Context context) {
		return basicReport.constructMapReport(datum, genome, chrName, context);
	}
	
	public void regionCoverReport(int depth, int noPCRdepth) {
		regionCoverReport.add(depth);
		rmdupRegionCoverReport.add(noPCRdepth);
	}
	
	public void finalizeUnmappedReport(String chrName) {
		if(options.isOutputUnmapped()) {
			unmappedReport.finalize(unmappedReport.getUnmappedSites(chrName));
		}
	}
	
	public void singleRegionReports(String chrName, long winStart, int winSize , PositionDepth pd) {
		if(options.getSingleRegion() != null)
			cnvSingleRegionReport.getStatisticString(chrName, (int) winStart, winSize, pd, "cnv");
	}
	
	public void insertReport(SamRecordDatum datum) {
		if((datum.getFlag() & 0x40) != 0) {
			int insert = datum.getInsertSize();
			if(Math.abs(insert) < 2000) {
				insertSize[Math.abs(insert)]++;
				if(!datum.isDup()) {
					insertSizeWithoutDup[Math.abs(insert)]++;
				}
			}
		}
	}
	
	protected void insertSizeReportReducerString(StringBuffer info, int[] insertSize) {
		for(int i = 0; i < insertSize.length; i++) {
			if(insertSize[i] != 0) {
				info.append(i);
				info.append("\t");
				info.append(insertSize[i]);
				info.append("\n");
			}
		}
	}
	
	public abstract void constructDepthReport(PositionDepth pd, int i, String chrName, long pos);

	public abstract String toReducerString(String sample, String chrName, boolean unmappedRegion);

	public void parseReport(LineReader lineReader, Text line, ReferenceShare genome) throws IOException {
		String lineString = line.toString();
		String chrName = "";
		if(lineString.contains("chrName:")) {
			String[] sampleSplit = line.toString().split(":");
			chrName = sampleSplit[1];
		}
		
		if(lineString.contains("Basic Information")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				basicReport.parse(line.toString());
			}
		}
		
		if(lineString.startsWith("cnv part single Region Statistic")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				cnvSingleRegionReport.parseReducerOutput(line.toString(), true);
			}
		}
		if(lineString.startsWith("cnv single Region Statistic")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				cnvSingleRegionReport.parseReducerOutput(line.toString(), false);
			}
		}
		if(lineString.startsWith("Region Depth")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				regionCoverReport.parseReducerOutput(line.toString());
			}
		}
		if(lineString.startsWith("RMDUP Region Depth")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				rmdupRegionCoverReport.parseReducerOutput(line.toString());
			}
		}
		
		if(lineString.contains("insert size information")) {
			fillInsertSize(lineReader, line, insertSize);
		}
		
		if(lineString.contains("insert size without dup information")) {
			fillInsertSize(lineReader, line, insertSizeWithoutDup);
		}
		
		if(lineString.contains("unmapped site information") && options.isOutputUnmapped()) {
			String[] splitArray = null;
			ArrayList<Long> unmappedSites =  unmappedReport.getUnmappedSites(chrName);
			while(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				if(line.toString().contains("unmapped site information")) {
					break;
				}
				splitArray = line.toString().split("\t");
				
				unmappedSites.add(Long.parseLong(splitArray[0]));
				unmappedSites.add(Long.parseLong(splitArray[1]));
			}
		}
	}
	
	public void write(FileSystem fs, String sampleName) throws IOException {
		if(cnvSingleRegionReport != null) {
			Map<Regiondata, SingleRegionStatistic> result = cnvSingleRegionReport.getResult();
			cnvSingleRegionReport.updateAllRegionAverageDeepth();

//			//输出小区域每个位点的深度，should setStatPosDepth(true) for cnvSingleRegionReport
//			String depthFileStr = options.getOutputPath() +
//					"/" + sampleName + ".region.depth.tsv.gz";
//			Path depthFilePath = new Path(depthFileStr);
//			FSDataOutputStream cnvDepthStream = fs.create(depthFilePath);
//			CompressionCodecFactory codecFactory = new CompressionCodecFactory(fs.getConf());
//			CompressionCodec codec = codecFactory.getCodec(depthFilePath);
//			CompressionOutputStream depthCompressedOutput = codec.createOutputStream(cnvDepthStream);
//			StringBuilder sb = new StringBuilder();
//			sb.append("#Chr\tPos\tRaw Depth\tRmdup depth\n");
//			for(Regiondata regionData : cnvSingleRegionReport.getRegion().getRegions()) {
//				SingleRegionStatistic singleRegionStat = result.get(regionData);
//				for (int i = 0; i < regionData.size(); i++) {
//					sb.append(regionData.getChrName());
//					sb.append("\t");
//					sb.append(i + regionData.getStart() + 1);
//					sb.append("\t");
//					sb.append(singleRegionStat.getPosDepth(i));
//					sb.append("\t");
//					sb.append(singleRegionStat.getPosRmdupDepth(i));
//					sb.append("\n");
//				}
//			}
//			depthCompressedOutput.write(sb.toString().getBytes());
//			depthCompressedOutput.close();
//			cnvDepthStream.close();


			String singleRegionFileStr = options.getOutputPath() +
					"/" +
					sampleName +
					".region.cov.tsv.gz";
			Path singleRegionPath = new Path(singleRegionFileStr);
			FSDataOutputStream singleRegionStream = fs.create(singleRegionPath);
			CompressionCodecFactory covCodecFactory = new CompressionCodecFactory(fs.getConf());
			CompressionCodec covCodec = covCodecFactory.getCodec(singleRegionPath);
			CompressionOutputStream covCompressedOutput = covCodec.createOutputStream(singleRegionStream);
			covCompressedOutput.write(SingleRegionStatistic.toReportTitleString().getBytes());
			covCompressedOutput.write('\n');

			String unsingleRegionFilePath = options.getOutputPath() +
					"/" +
					sampleName +
					".region.lowdepth.cov.tsv";
			Path unsingleRegionPath = new Path(unsingleRegionFilePath);
			FSDataOutputStream unsingleRegionwriter = fs.create(unsingleRegionPath);
			unsingleRegionwriter.write(SingleRegionStatistic.toReportTitleString().getBytes());
			unsingleRegionwriter.write('\n');

			for(Regiondata regionData : cnvSingleRegionReport.getRegion().getRegions()) {
				SingleRegionStatistic singleRegionStat = result.get(regionData);
				if(singleRegionStat.getDepth(regionData) > options.getMinSingleRegionDepth())
					covCompressedOutput.write(result.get(regionData).toReportString(regionData, cnvSingleRegionReport.getAllRegionAverageDepth()).getBytes());
				else {
					covCompressedOutput.write(result.get(regionData).toReportString(regionData, cnvSingleRegionReport.getAllRegionAverageDepth()).getBytes());
					unsingleRegionwriter.write(result.get(regionData).toReportString(regionData, cnvSingleRegionReport.getAllRegionAverageDepth()).getBytes());
				}
			}
			covCompressedOutput.close();
			unsingleRegionwriter.close();
		}

		if(regionCoverReport != null) {
			String RegionDepthFilePath = options.getOutputPath() +	"/" + sampleName + ".depth.txt";
			Path regionDepthPath = new Path(RegionDepthFilePath);
			FSDataOutputStream regionDepthWriter = fs.create(regionDepthPath);
			int[] depth = regionCoverReport.getDepthArray();
			
			if(rmdupRegionCoverReport != null)
			{
				int[] rmdepth = rmdupRegionCoverReport.getDepthArray();

				StringBuilder sb = new StringBuilder();
				for(int i = 0; i < depth.length; i++) {
			             sb.append(i);
			             sb.append("\t");
				     sb.append(depth[i]);
			             sb.append("\t");
				     sb.append(rmdepth[i]);
			             sb.append("\n");
				}
				regionDepthWriter.write(sb.toString().getBytes());
			}else{
				regionDepthWriter.write(regionCoverReport.toString().getBytes());
			}
			regionDepthWriter.close();
		}
		
		StringBuilder reportFilePath = new StringBuilder();
//		reportFilePath.append(options.getOutputPath());
//		reportFilePath.append("/");
//		reportFilePath.append(sampleName);
//		reportFilePath.append(".unmapped.bed");
//		Path bedPath = new Path(reportFilePath.toString());
//		FSDataOutputStream bedwriter = fs.create(bedPath);
//
//		StringBuilder bedString = new StringBuilder();
//		for(String chrName:unmappedReport.getUnmappedSites().keySet()) {
//			ArrayList<Long> sites = unmappedReport.getUnmappedSites(chrName);
//			for(int i = 0; i < sites.size(); i += 2) {
//				bedString.append(chrName);
//				bedString.append("\t");
//				bedString.append(sites.get(i));
//				bedString.append("\t");
//				bedString.append(sites.get(i+1));
//				bedString.append("\n");
//			}
//		}
//		bedwriter.write(bedString.toString().getBytes());
//		bedwriter.close();
		
		reportFilePath.setLength(0);
		reportFilePath.append(options.getOutputPath());
		reportFilePath.append("/");
		reportFilePath.append(sampleName);
		reportFilePath.append(".insert.xls");
		Path insertPath = new Path(reportFilePath.toString());
		FSDataOutputStream insertwriter = fs.create(insertPath);
		
		StringBuilder insertString = new StringBuilder();
		for(int i = 0; i < insertSize.length; i++) {
			insertString.append(i);
			insertString.append("\t");
			insertString.append(insertSize[i]);
			insertString.append("\n");
		}
		insertwriter.write(insertString.toString().getBytes());
		insertwriter.close();
		
		reportFilePath.setLength(0);
		reportFilePath.append(options.getOutputPath());
		reportFilePath.append("/");
		reportFilePath.append(sampleName);
		reportFilePath.append(".insert_without_dup.xls");
		Path insertWithoutDupPath = new Path(reportFilePath.toString());
		FSDataOutputStream insertWithoutDupwriter = fs.create(insertWithoutDupPath);
		
		StringBuilder insertWithoutDupString = new StringBuilder();
		for(int i = 0; i < insertSizeWithoutDup.length; i++) {
			insertWithoutDupString.append(i);
			insertWithoutDupString.append("\t");
			insertWithoutDupString.append(insertSizeWithoutDup[i]);
			insertWithoutDupString.append("\n");
		}
		insertWithoutDupwriter.write(insertWithoutDupString.toString().getBytes());
		insertWithoutDupwriter.close();
	}
	
	private void fillInsertSize(LineReader lineReader, Text line, int[] insertSize) throws RuntimeException, IOException {
		String[] splitArray = null;
		while(lineReader.readLine(line) > 0 && line.getLength() != 0) {
			if(line.toString().contains("insert size")) {
				break;
			}
			splitArray = line.toString().split("\t");
			int index = Integer.parseInt(splitArray[0]);
			insertSize[index] += Integer.parseInt(splitArray[1]);
		}
	}
	
	public SingleRegionReport getCNVSingleRegionReport() {
		return cnvSingleRegionReport;
	}
	
	public ReferenceShare getReference() {
		return genome;
	}
}
