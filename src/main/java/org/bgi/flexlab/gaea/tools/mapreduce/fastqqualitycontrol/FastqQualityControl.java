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

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.lib.input.MultipleInputs;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.MultipleOutputs;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.adaptor.AdaptorInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.fastq.FastqInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.fastq.FastqMultipleSample;
import org.bgi.flexlab.gaea.data.mapreduce.input.fastq.FastqRecordReader;
import org.bgi.flexlab.gaea.data.mapreduce.input.fastq.FastqSample;
import org.bgi.flexlab.gaea.data.structure.reads.report.FastqQualityControlReporterIO;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.PairEndAggregatorMapper;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;

import java.util.Map;

public class FastqQualityControl extends ToolsRunner {

	public FastqQualityControl() {
		this.toolsDescription = "Gaea fastq quality control\n";
	}

	@Override
	public int run(String[] args) throws Exception {
		BioJob job = BioJob.getInstance();
		Configuration conf = job.getConfiguration();
		String[] remainArgs =  remainArgs(args,conf);
		
		FastqQualityControlOptions option = new FastqQualityControlOptions();
		option.parse(remainArgs);
		conf.setInt(FastqRecordReader.READ_NAME_TYPE, option.getReadType());
		option.setHadoopConf(remainArgs, conf);

		job.setJobName("GaeaFastqQC");
		job.setJarByClass(FastqQualityControl.class);
		job.setMapperClass(PairEndAggregatorMapper.class);
		job.setReducerClass(FastqQualityControlReducer.class);

		job.setInputFormatClass(FastqInputFormat.class);
		job.setOutputFormatClass(TextOutputFormat.class);
		job.setNumReduceTasks(option.getReducerNumber());
		job.setOutputKeyValue(Text.class, Text.class, NullWritable.class,
				Text.class);

		FastqMultipleSample sample = null;
		if (option.getMultiSampleList() != null
				&& option.getMultiSampleList() != "") {
			sample = new FastqMultipleSample(option.getMultiSampleList(), true);
			Map<String, FastqSample> sampleList = sample.getSampleList();

			for (FastqSample sl : sampleList.values()) {
				if (sl.getFastq1() != null) {
					MultipleInputs.addInputPath(job, new Path(sl.getFastq1()),
							FastqInputFormat.class);
				} else {
					System.err.println(sl.getSampleName() + " has no fq1!");
					System.exit(1);
				}
				if (sl.getFastq2() != null) {
					MultipleInputs.addInputPath(job, new Path(sl.getFastq2()),
							FastqInputFormat.class);
				} else {
					System.err.println(sl.getSampleName() + " is SE data!");
				}
				if (sl.getAdapter1() != null) {
					MultipleInputs.addInputPath(job,
							new Path(sl.getAdapter1()),
							AdaptorInputFormat.class);
				}
				if (sl.getAdapter2() != null) {
					MultipleInputs.addInputPath(job,
							new Path(sl.getAdapter2()),
							AdaptorInputFormat.class);
				}
			}
		} else {
			if (option.getInputFastq1() != null) {
				MultipleInputs.addInputPath(job,
						new Path(option.getInputFastq1()),
						FastqInputFormat.class);
			}
			if (option.getInputFastq2() != null) {
				MultipleInputs.addInputPath(job,
						new Path(option.getInputFastq2()),
						FastqInputFormat.class);
			}
			if (option.getAdapter1() != null) {
				MultipleInputs.addInputPath(job,
						new Path(option.getAdapter1()),
						AdaptorInputFormat.class);
			}
			if (option.getAdapter2() != null) {
				MultipleInputs.addInputPath(job,
						new Path(option.getAdapter2()),
						AdaptorInputFormat.class);
			}
		}

		Path outputPath = new Path(option.getOutputDirectory() + "/out_fq");
		FileOutputFormat.setOutputPath(job, outputPath);
		MultipleOutputs.addNamedOutput(job, "filterStatistic",
				TextOutputFormat.class, NullWritable.class, Text.class);
		MultipleOutputs.addNamedOutput(job, "qualFreqStatistic",
				TextOutputFormat.class, NullWritable.class, Text.class);

		if (job.waitForCompletion(true)) {
			FastqQualityControlReporterIO report = new FastqQualityControlReporterIO(
					sample, option.isMultiStatis());
			report.mergeReport(outputPath, conf,
					new Path(option.getOutputDirectory()));
			return 0;
		} else {
			return 1;
		}
	}
}