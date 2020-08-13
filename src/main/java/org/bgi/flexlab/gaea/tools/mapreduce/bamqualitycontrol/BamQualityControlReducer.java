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

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.data.structure.positioninformation.CompoundInformation;
import org.bgi.flexlab.gaea.data.structure.positioninformation.depth.PositionDepth;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.RegionResultReport;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.ReportBuilder;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.ResultReport;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.WholeGenomeResultReport;
import org.bgi.flexlab.gaea.util.SamRecordDatum;

import java.io.IOException;

public class BamQualityControlReducer extends Reducer<Text, Text, NullWritable, Text>{
	
	private BamQualityControlOptions options;
		
	private ResultReport reportType;
	
	private ReportBuilder reportBuilder;
		
	private PositionDepth deep;
				
	@Override
	protected void setup(Context context) throws IOException {
		options = new BamQualityControlOptions();
		Configuration conf = context.getConfiguration();
		options.getOptionsFromHadoopConf(conf);
		
		reportBuilder = new ReportBuilder();
		if ((options.getRegion() != null) || (options.getBedfile() != null))
			reportType = new RegionResultReport(options, conf);
		else
			reportType = new WholeGenomeResultReport(options);
	
	}
	
	@Override
	public void reduce(Text key, Iterable<Text> values,Context context) throws IOException, InterruptedException {
		String[] keySplit = key.toString().split(":");	
		String sampleName = keySplit[0];	
		String chrName = keySplit[1];					
		long winNum = Long.parseLong(keySplit[2]);		
		ChromosomeInformationShare chrInfo = null;
		try{
			chrInfo = reportType.getReference().getChromosomeInfo(chrName);
		} catch (Exception e) {
			if(chrName.equals("-1")) {
//				do nothing, it means unmapped area
			} else {
				throw new RuntimeException(e);
			}
		}
		reportBuilder.setReportChoice(reportType);
		reportBuilder.initReports(sampleName, chrName);

		if(reportBuilder.unmappedReport(winNum, chrName, values)) {
			ResultReport report = reportBuilder.build();
			context.write(NullWritable.get(), new Text(report.toReducerString(sampleName, chrName, true)));
			return;
		}

		long start = winNum * BamQualityControl.WINDOW_SIZE; // 窗口起始坐标

		if (chrInfo == null)
			return;

		if(start > chrInfo.getLength()) {
			context.getCounter("Exception", "window start beyond chromosome").increment(1);
			return;
		}

		long readPos;
		long winEnd = start + BamQualityControl.WINDOW_SIZE - 1;
		
		int winSize = BamQualityControl.WINDOW_SIZE;
		if(winEnd > chrInfo.getLength()) {
			winSize = (int) (chrInfo.getLength() - start);
			winEnd = chrInfo.getLength() - 1;
			context.getCounter("Exception", "window end > chr length:" ).increment(1);
		}
		//position depth
		deep = new PositionDepth(winSize, options.isGenderDepth(), reportBuilder.getSampleLaneSzie(sampleName));
				
		for(Text value : values) {
			SamRecordDatum datum = new SamRecordDatum();

			if(!datum.parseBamQC(value.toString())) {
				context.getCounter("Exception", "parse mapper output error").increment(1);
				continue;
			}
			readPos = datum.getPosition(); 
			if (readPos < 0) {
				context.getCounter("Exception", "read start pos less than zero").increment(1);
				continue;
			}
			
			if(!deep.add(new CompoundInformation<SamRecordDatum>((int) start, winSize, datum, chrInfo)))
				context.getCounter("Exception", "null read info in depth class").increment(1);
			
			if(datum.isRepeat())
				continue;
			
			if (!reportBuilder.mappedReport(datum, chrName, context)) 
				continue;
				
			reportBuilder.insertReport(datum);
		}
				
		for(int i = 0; i < winSize; i++) {
			long pos = start + i;
			int depth = deep.getPosDepth(i);
			if(depth < 0) {
				throw new RuntimeException("sample:" + sampleName + " chr:" + chrName + " windSize:" + winSize + " position:" + i + " winStart:" + start + " depth:" + depth + " < 0.");
			}
			reportBuilder.constructDepthReport(deep, i, chrName, pos);
		}
		
		//last unmapped sites
		//BED format: 0,1 stand for 1 position
		if(options.isOutputUnmapped()) {
			reportBuilder.finalizeUnmappedReport(chrName);
		}
		//sum single region
		reportBuilder.singleRegionReports(chrName, start, winSize, deep);
		
		//write reducer
		ResultReport report = reportBuilder.build();
		context.write(NullWritable.get(), new Text(report.toReducerString(sampleName, chrName, false)));
	} 

}
