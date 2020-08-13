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
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.bam.GaeaAnySAMInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.cram.GaeaCramInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.BamReport;

import static org.bgi.flexlab.gaea.data.mapreduce.input.cram.GaeaCramRecordReader.INPUTFORMAT_REFERENCE;

public class BamQualityControl extends ToolsRunner{
	
	public final static int WINDOW_SIZE = 200000;
	
	public BamQualityControl() {
			this.toolsDescription = "Gaea bam quality control\n"
					+ "The purpose of bam quality control is to attain statistics information"
					+ "of the bam file";
	}
	 
	private BamQualityControlOptions options;
	
	@Override
	public int run(String[] args) throws Exception {
		BioJob job = BioJob.getInstance();
		Configuration conf = job.getConfiguration();
		String[] remainArgs =  remainArgs(args,conf);
		options = new BamQualityControlOptions();
		options.parse(remainArgs);
		options.setHadoopConf(remainArgs, conf);
		
		if(options.isDistributeCache()) {
			ReferenceShare.distributeCache(options.getReferenceSequencePath(), job);
		}

		boolean iscram = options.getInputs().get(0).toString().endsWith("cram");
		if(iscram){
			if(options.getLocalReferenceSequencePath() == null)
				throw new RuntimeException("Please set local reference for cram (-f).");
			conf.set(INPUTFORMAT_REFERENCE, options.getLocalReferenceSequencePath());
		}

		SamHdfsFileHeader.loadHeader(options.getInputs(), conf, new Path(options.getOutputPath()), iscram);

		job.setJobName("BamQualityControl");
		job.setJarByClass(BamQualityControl.class);
		job.setMapperClass(BamQualityControlMapper.class);
		job.setReducerClass(BamQualityControlReducer.class);
		job.setOutputKeyValue(Text.class, Text.class, 
				NullWritable.class, Text.class);
		job.setNumReduceTasks(options.getReducerNum());

		FileInputFormat.setInputPaths(job, options.getInputs().toArray(new Path[options.getInputs().size()]));
		if(iscram)
			job.setInputFormatClass(GaeaCramInputFormat.class);
		else
			job.setInputFormatClass(GaeaAnySAMInputFormat.class);

		job.setOutputFormatClass(TextOutputFormat.class);
		
		FileOutputFormat.setOutputPath(job, new Path(options.getTempPath()));
		if(job.waitForCompletion(true)) {
			BamReport.getOutput(options, conf, new Path(options.getTempPath()));
			return 0;
		} else {
			return 1;
		}
	}
	
	public static void main(String[] args) throws Exception {
		BamQualityControl bamqc = new BamQualityControl();
		bamqc.run(args);
	}
}
