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

import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.structure.positioninformation.depth.PositionDepth;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty;
import org.bgi.flexlab.gaea.tools.mapreduce.bamqualitycontrol.BamQualityControlOptions;

import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;

public class WholeGenomeResultReport extends ResultReport{
			
	private Map<String, WholeGenomeCoverReport> coverReports;
		
	public WholeGenomeResultReport(BamQualityControlOptions options) throws IOException {
		super(options);
	}
	
	public void initReports() throws IOException {
		super.initReports();
		coverReports = new ConcurrentHashMap<String, WholeGenomeCoverReport>();
	}
	
	public void initReports(String chrName) throws IOException {
		super.initReports();
		if(chrName.equals("-1"))
			return;
		else {
			coverReports = new ConcurrentHashMap<>();
			ChromosomeInformationShare chrInfo = genome.getChromosomeInfo(chrName);
			coverReports.put(chrName, new WholeGenomeCoverReport(chrInfo));
		}
	}
	
	@Override
	public void constructDepthReport(PositionDepth pd, int i, String chrName, long pos) {
		int depth = pd.getPosDepth(i);
		int noPCRdepth = pd.getRMDupPosDepth(i);
		super.regionCoverReport(depth, noPCRdepth);
		if(options.isOutputUnmapped() && depth != 0 )
			unmappedReport.updateUnmappedSites(pos, unmappedReport.getUnmappedSites(chrName));
		coverReports.get(chrName).constructDepthReport(pd, i, (int)pos);
	}
	
	@Override
	public String toReducerString(String sample, String chrName, boolean unmappedRegion) {
		StringBuffer info = new StringBuffer();
		
		if(chrName.equals("-1")) {
			info.append("sample:");
			info.append(sample);
			info.append("\n");
			info.append(basicReport.toReducerString());
			return info.toString();
		}
		info.append("sample:");
		info.append(sample);
		info.append("\n");
		info.append("chrName:");
		info.append(chrName);
		info.append("\n");
		info.append(basicReport.toReducerString());
		if(!unmappedRegion) {
			info.append(regionCoverReport.toReducerString());
			info.append(rmdupRegionCoverReport.toReducerString("RMDUP Region Depth"));
			for(String key : coverReports.keySet()) {
				WholeGenomeCoverReport cover = coverReports.get(key);
				info.append(cover.toReducerString(key));
			}
			info.append("insert size information:\n");
			insertSizeReportReducerString(info, insertSize);
			info.append("insert size information\n");
			
			info.append("insert size without dup information:\n");
			insertSizeReportReducerString(info, insertSizeWithoutDup);
			info.append("insert size without dup information\n");
			
			info.append("unmapped site information:\n");
			unmappedReport.toReducerString();
			info.append("unmapped site information\n");
			if(cnvSingleRegionReport != null)
				info.append(cnvSingleRegionReport.toReducerString());
		}
		return info.toString();
	}
	
	@Override
	public void parseReport(LineReader lineReader, Text line, ReferenceShare genome) throws IOException {
		super.parseReport(lineReader, line, genome);
		String lineString = line.toString();
		if(lineString.contains("Cover Information")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				String[] splitArray = line.toString().split("\t");
				WholeGenomeCoverReport coverReport = null;
				for(String keyValue : splitArray) {
					if(keyValue.split(" ").length == 1) {
						String chrName = keyValue;
						if(!coverReports.containsKey(chrName)) {
							ChromosomeInformationShare chrInfo = genome.getChromosomeInfo(chrName);
							coverReport = new WholeGenomeCoverReport(chrInfo);
							coverReports.put(chrName, coverReport);
						} else {
							coverReport = coverReports.get(chrName);
						}
					} else {
						assert coverReport != null;
						coverReport.parse(keyValue, genome);
					}
				}
			}
		}
	}
	
	@Override
	public void write(FileSystem fs, String sampleName) throws IOException {
		DecimalFormat df = new DecimalFormat("0.000");
		df.setRoundingMode(RoundingMode.HALF_UP);

		super.write(fs, sampleName);
		StringBuffer reportFilePath = new StringBuffer();
		reportFilePath.append(options.getOutputPath());
		reportFilePath.append("/");
		reportFilePath.append(sampleName);
		reportFilePath.append(".bam.report.txt");
		Path reportPath = new Path(reportFilePath.toString());
		FSDataOutputStream reportwriter = fs.create(reportPath);
		StringBuffer info = new StringBuffer();
		info.append(basicReport.toString());

		long depth = 0;
		long rmdupDepth = 0;
		long coverBaseNum = 0;
		long fourxCoverBaseNum = 0;
		long tenxCoverBaseNum = 0;
		long tenxNonnCoverBaseNum = 0;
		long thirtyxCoverBaseNum = 0;
		long thirtyxNonnCoverBaseNum = 0;
		long nonnCoverBaseNum = 0;
		long fourxNonnCoverBaseNum = 0;
		long genomeLength = 0;
		long nonnGenomeLength = 0;
		TreeSet<String> keys = new TreeSet<String>(coverReports.keySet());
		for(String key : keys) {
			WholeGenomeCoverReport cover = coverReports.get(key);
			coverBaseNum += cover.getCoverBaseNum();
			fourxCoverBaseNum += cover.getFourXCoverBaseNum();
			fourxNonnCoverBaseNum += cover.getFourXNonNCoverBaseNum();
			tenxCoverBaseNum += cover.getTenXCoverBaseNum();
			tenxNonnCoverBaseNum += cover.getTenXNonNCoverBaseNum();
			thirtyxCoverBaseNum += cover.getThirtyXCoverBaseNum();
			thirtyxNonnCoverBaseNum += cover.getThirtyXNonnCoverBaseNum();
			nonnCoverBaseNum += cover.getNonNCoverBaseNum();
			genomeLength += cover.getChrInfo().getLength();
			nonnGenomeLength += cover.getChrInfo().getNonNbaselength();
			depth += cover.getDepth();
			rmdupDepth += cover.getRmdupDepth();
		}

		info.append("Average depth:\t");
		info.append(df.format(depth/(double)coverBaseNum));
		info.append("\nAverage depth(rmdup):\t");
		info.append(df.format(rmdupDepth/(double)coverBaseNum));
		info.append("\nCoverage (>=1x):\t");
		info.append(df.format(100*coverBaseNum/(double)genomeLength));
		info.append("%\nCoverage (>=4x):\t");
		info.append(df.format(100*fourxCoverBaseNum/(double)genomeLength));
		info.append("%\nCoverage (>=10x):\t");
		info.append(df.format(100*tenxCoverBaseNum/(double)genomeLength));
		info.append("%\nCoverage (>=30x):\t");
		info.append(df.format(100*thirtyxCoverBaseNum/(double)genomeLength));
		info.append("%\nNonN Average depth:\t");
		info.append(df.format(depth/(double)nonnCoverBaseNum));
		info.append("\nNonN Average depth(rmdup):\t");
		info.append(df.format(rmdupDepth/(double)nonnCoverBaseNum));
		info.append("\nNonN Coverage (>=1x):\t");
		info.append(df.format(100*nonnCoverBaseNum/(double)nonnGenomeLength));
		info.append("%\nNonN Coverage (>=4x):\t");
		info.append(df.format(100*fourxNonnCoverBaseNum/(double)nonnGenomeLength));
		info.append("%\nNonN Coverage (>=10x):\t");
		info.append(df.format(100*tenxNonnCoverBaseNum/(double)nonnGenomeLength));
		info.append("%\nNonN Coverage (>=30x):\t");
		info.append(df.format(100*thirtyxNonnCoverBaseNum/(double)nonnGenomeLength));
		if(cnvSingleRegionReport != null) {
			info.append("%\n[Target] Len of region(without XY):\t");
			info.append((long) cnvSingleRegionReport.getRegionSizeTotal());
			info.append("\n[Target] Average depth:\t");
			info.append(df.format(cnvSingleRegionReport.getAllRegionAverageDepth()));
			info.append("\n[Target] Average depth(rmdup):\t");
			info.append(df.format(cnvSingleRegionReport.getAllRegionAverageRmdupDepth()));
			info.append("\n[Target] Coverage (>=1x):\t");
			info.append(df.format(100*cnvSingleRegionReport.getAllRegionCoverage()));
			info.append("%\n[Target] Coverage (>=4x):\t");
			info.append(df.format(100*cnvSingleRegionReport.getAllRegionFourxCoverage()));
			info.append("%\n[Target] Coverage (>=10x):\t");
			info.append(df.format(100*cnvSingleRegionReport.getAllRegionTenxCoverage()));
			info.append("%\n[Target] Coverage (>=30x):\t");
			info.append(df.format(100*cnvSingleRegionReport.getAllRegionThirtyxCoverage()));
		}
		info.append("%\n\ncoverage information:\n");
		TreeSet<String> keys2 = new TreeSet<String>(coverReports.keySet());
		for(String key : keys2) {
			WholeGenomeCoverReport cover = coverReports.get(key);
			info.append(cover.toString());
		}
		reportwriter.write(info.toString().getBytes());
		reportwriter.close();
	}
}
