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

import htsjdk.samtools.SAMFileHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.MultipleOutputs;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.bgi.flexlab.gaea.data.exception.FileNotExistException;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.output.bam.GaeaBamOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.util.HdfsFileManager;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.WindowsBasedSamRecordMapper;
import org.bgi.flexlab.gaea.tools.recalibrator.report.RecalibratorReportTableEngine;
import org.bgi.flexlab.gaea.tools.recalibrator.table.RecalibratorTableCombiner.NonRecalibratorPathFilter;
import org.seqdoop.hadoop_bam.SAMFormat;

import java.io.IOException;
import java.util.ArrayList;

public class Realigner extends ToolsRunner {

	public Realigner() {
		this.toolsDescription = "Gaea realigner and base quality recalibrator!\n";
	}

	public final static String RECALIBRATOR_REPORT_TABLE_NAME = "bqsr.report.table";

	private final static SAMFormat format = SAMFormat.BAM;

	private RealignerExtendOptions options = null;
	private RealignerOptions option = null;
	private SAMFileHeader header = null;

	private int runRealigner(String[] args) throws IOException, ClassNotFoundException, InterruptedException {
		BioJob job = BioJob.getInstance();
		Configuration conf = job.getConfiguration();
		String[] remainArgs = remainArgs(args, conf);

		options = new RealignerExtendOptions();
		options.parse(remainArgs);

		option = options.getRealignerOptions();
		
		String jobName = "Gaea realigner and recalibrator";

		if (options.isRecalibration() && !options.isRealignment()) {
			job.setOnlyBaseRecalibrator(true);
			jobName = "GaeaRecalibrator";
		}else if(options.isRealignment() && !options.isRecalibration())
			jobName = "GaeaRealigner";
		
		if(option.isMultiSample())
			job.setMultipleSample();

		job.setJobName(jobName);
		
		option.setHadoopConf(remainArgs, conf);

		header = job.setHeader(new Path(option.getRealignerInput()), new Path(options.getCommonOutput()));

		job.setAnySamInputFormat(option.getInputFormat());
		job.setOutputFormatClass(GaeaBamOutputFormat.class);
		job.setOutputKeyValue(WindowsBasedWritable.class, SamRecordWritable.class, NullWritable.class,
				SamRecordWritable.class);

		job.setJarByClass(Realigner.class);
		job.setWindowsBasicMapperClass(WindowsBasedSamRecordMapper.class, option.getWindowsSize(),option.getExtendSize());
		job.setReducerClass(RealignerReducer.class);
		job.setNumReduceTasks(option.getReducerNumber());

		FileInputFormat.setInputPaths(job, new Path(option.getRealignerInput()));
		FileOutputFormat.setOutputPath(job, new Path(option.getRealignerOutput()));

		if (options.isRecalibration())
			MultipleOutputs.addNamedOutput(job, RecalibratorContextWriter.RECALIBRATOR_TABLE_TAG,
					TextOutputFormat.class, NullWritable.class, Text.class);

		if (job.waitForCompletion(true)) {
			if (options.isRecalibration())
				return mergeReportTable(options.getBqsrOptions(), header,
						options.getCommonOutput() + RECALIBRATOR_REPORT_TABLE_NAME);
			return 0;
		}

		return 1;
	}

	private int runFixMate(String[] args) throws IOException, ClassNotFoundException, InterruptedException {
		BioJob job = BioJob.getInstance();
		
		String jobName = "Gaea fixmate and print reads";

		Configuration conf = job.getConfiguration();

		if (options == null) {
			String[] remainArgs = remainArgs(args, conf);

			options = new RealignerExtendOptions();
			options.parse(remainArgs);

			option = options.getRealignerOptions();
		}

		// set bqsr table path
		if (options.isRecalibration())
			conf.set(RECALIBRATOR_REPORT_TABLE_NAME, options.getCommonOutput() + RECALIBRATOR_REPORT_TABLE_NAME);

		String[] remainArgs = remainArgs(args, conf);
		option.setHadoopConf(remainArgs, conf);

		job.setHeader(options.getCommonOutput() + "/" + SamHdfsFileHeader.BAM_HEADER_FILE_NAME);

		job.setAnySamInputFormat(format);
		job.setOutputFormatClass(GaeaBamOutputFormat.class);

		job.setJarByClass(Realigner.class);
		job.setMapperClass(FixmateMapper.class);
		job.setReducerClass(FixmateReducer.class);

		if (!options.isRealignment()) {
			jobName = "GaeaPrintReads";
			job.setNumReduceTasks(0);
			job.setOutputKeyValue(NullWritable.class, SamRecordWritable.class, NullWritable.class,
					SamRecordWritable.class);
		} else {
			job.setOutputKeyValue(Text.class, SamRecordWritable.class, NullWritable.class, SamRecordWritable.class);
			job.setNumReduceTasks(option.getReducerNumber());
			
			if(!options.isRecalibration())
				jobName = "GaeaFixmate";
		}
		
		job.setJobName(jobName);
		
		ArrayList<Path> lists = inputFilter(option.getFixmateInput(options.isRealignment()),conf,new NonRecalibratorPathFilter());
		
		for(Path p : lists){
			FileInputFormat.addInputPath(job, p);
		}
		lists.clear();

		FileOutputFormat.setOutputPath(job, new Path(option.getFixmateOutput()));

		return job.waitForCompletion(true) ? 0 : 1;
	}
	
	private ArrayList<Path> inputFilter(String path,Configuration conf,NonRecalibratorPathFilter filter){
		ArrayList<Path> lists = new ArrayList<Path>();
		
		Path p = new Path(path);
		
		FileSystem fs = HdfsFileManager.getFileSystem(p, conf);
		FileStatus status = null;
		try {
			status = fs.getFileStatus(p);
		} catch (IOException e2) {
			throw new FileNotExistException(p.getName());
		}
		
		if(status.isFile()){
			if(filter.accept(status.getPath()))
				lists.add(status.getPath());
		}else{
			FileStatus[] stats = null;
			try{
				stats = fs.listStatus(p,filter);
			}catch(IOException e){
				throw new RuntimeException(e.toString());
			}

			for (FileStatus file : stats) {
				if(!file.isFile()){
					throw new RuntimeException("input directory cann't contains sub directory!!");
				}else{
					if(filter.accept(file.getPath()))
						lists.add(file.getPath());
				}
			}
		}
		
		if(lists.size() == 0)
			return null;
		return lists;
	}

	private int mergeReportTable(RecalibratorOptions bqsrOption, SAMFileHeader header, String output) {
		RecalibratorHdfsReportWriter writer = new RecalibratorHdfsReportWriter(output);
		RecalibratorReportTableEngine engine = new RecalibratorReportTableEngine(bqsrOption, header, writer);
		engine.writeReportTable(option.getRealignerOutput());
		return 0;
	}

	@Override
	public int run(String[] args) throws Exception {
		int res = runRealigner(args);

		if (res != 0) {
			throw new RuntimeException("Realigner is failed!");
		}

		res = runFixMate(args);

		return res;
	}
}
