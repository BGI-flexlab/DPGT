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
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer.Context;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.structure.positioninformation.depth.PositionDepth;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.util.SamRecordDatum;

import java.io.IOException;

public class ReportBuilder {
	
	ResultReport report;
	
	public void setReportChoice(ResultReport report) {
		this.report = report;
	}
	
	public void initReports(String sampleName, String chrName) throws IOException {
		if(report instanceof WholeGenomeResultReport)
			((WholeGenomeResultReport) report).initReports(chrName);
		else if(report instanceof RegionResultReport)
			((RegionResultReport) report).initReports(sampleName);
	}
	
	public void initReports(String sampleName) throws IOException {
		if(report instanceof WholeGenomeResultReport)
			((WholeGenomeResultReport) report).initReports();
		else if(report instanceof RegionResultReport)
			((RegionResultReport) report).initReports(sampleName);
	}
	
	public boolean unmappedReport(long winNum, String chrName, Iterable<Text> values) {
		return report.unmappedReport(winNum, chrName, values);
	}

	public boolean unmappedReport(long winNum, String chrName, SAMRecordIterator values) {
		return report.unmappedReport(winNum, chrName, values);
	}
	
	public void finalizeUnmappedReport(String chrName) {
		report.finalizeUnmappedReport(chrName);
	}
	
	public boolean mappedReport(SamRecordDatum datum, String chrName, Context context) {
		return report.mappedReport(datum, chrName, context);
	}
	
	public void constructDepthReport(PositionDepth pd, int i, String chrName, long pos) {
		report.constructDepthReport(pd, i, chrName, pos);
	}
	
	public void singleRegionReports(String chrName, long winStart, int winSize , PositionDepth pd) {
		report.singleRegionReports(chrName, winStart, winSize, pd);
	}
	
	public void insertReport(SamRecordDatum datum) {
		report.insertReport(datum);
	}
	
	public int getSampleLaneSzie(String sample) {
		return report instanceof RegionResultReport ? ((RegionResultReport) report).getSampleLaneSize(sample) : 0;
	}
	
	public void parseReport(LineReader lineReader, Text line, ReferenceShare genome) throws IOException {
		report.parseReport(lineReader, line, genome);
	}
	
	public ResultReport build() {
		return report;
	}
}
