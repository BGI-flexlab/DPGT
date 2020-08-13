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

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.tools.mapreduce.bamqualitycontrol.BamQualityControlOptions;

import java.io.IOException;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class BamReport {
	
	public static void getOutput(BamQualityControlOptions options, Configuration conf, Path oPath) throws IOException {
		ReportBuilder reportBuilder = new ReportBuilder();
		ResultReport reportType;
		ReferenceShare genome = new ReferenceShare();
		genome.loadChromosomeList(options.getReferenceSequencePath());
		Map <String, ResultReport> reports = new ConcurrentHashMap<>();

		FileSystem fs = oPath.getFileSystem(conf);
		FileStatus filelist[] = fs.listStatus(oPath);
		for(int i = 0; i < filelist.length; i++) {
			if(!filelist[i].isDirectory() && !filelist[i].getPath().toString().startsWith("_")) {
				FSDataInputStream reader = fs.open(filelist[i].getPath());
				LineReader lineReader = new LineReader(reader, conf);
				Text line = new Text();
				while(lineReader.readLine(line) > 0) {
					String lineString = line.toString();
					if(line.getLength() == 0)
						continue;

					if(lineString.contains("sample:")) {
						String sample = line.toString().split(":")[1];
						if(!reports.containsKey(sample)) {
							if ((options.getRegion() != null) || (options.getBedfile() != null))
								reportType = new RegionResultReport(options, conf);
							else
								reportType = new WholeGenomeResultReport(options);
							reports.put(sample, reportType);
							reportBuilder.setReportChoice(reportType);
							reportBuilder.initReports(sample);
						} else {
							reportType = reports.get(sample);
							reportBuilder.setReportChoice(reportType);
						}
					}
					reportBuilder.parseReport(lineReader, line, genome);					

				}
				lineReader.close();
				reader.close();
			}	
		}
		
		for(String sampleName:reports.keySet()) {
			System.err.println("sample:" + sampleName);
			ResultReport report = reports.get(sampleName);
			report.write(fs, sampleName);
		}
		fs.delete(oPath, true);
		fs.close();
	}
}
