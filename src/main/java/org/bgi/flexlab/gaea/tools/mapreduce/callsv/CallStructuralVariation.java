package org.bgi.flexlab.gaea.tools.mapreduce.callsv;

import java.io.IOException;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.bam.GaeaAnySAMInputFormat;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.tools.callsv.NewMapKey;
import org.bgi.flexlab.gaea.tools.callsv.ReduceGroupingComparator;
import org.bgi.flexlab.gaea.tools.callsv.SamWritable;

public class CallStructuralVariation extends ToolsRunner{
	
	/*
	public static void main(String[] args) throws ClassNotFoundException, IOException, InterruptedException {
		CallStructuralVariation sv = new CallStructuralVariation();
		int res = sv.runCallStructuralVariation(args);
		System.exit(res);
	}
	*/
	
	public CallStructuralVariation() {
		this.toolsDescription = "Gaea structural variantion calling";
	}

	private CallStructuralVariationOptions options = null;
	
	private int runCallStructuralVariation(String[] args) throws IOException, ClassNotFoundException, InterruptedException {
		
		/**
		 * set job1 info
		 */
		BioJob job = BioJob.getInstance();
		
		Configuration conf = job.getConfiguration();
		
		conf.set("mapreduce.reduce.shuffle.input.buffer.percent", "0.5");
		conf.set("mapreduce.reduce.shuffle.memory.limit.percent", "0.1f");
		conf.set("mapreduce.reduce.shuffle.parallelcopies", "10");
		
		String[] remainArgs1 = remainArgs(args, conf);
		options = new CallStructuralVariationOptions();
		options.parse(remainArgs1);
		options.setHadoopConf(remainArgs1, conf);
		
		job.setJobName("CallSV");
		job.setJarByClass(CallStructuralVariation.class);
		job.setMapperClass(CallStructuralVariationMapper.class);
		job.setReducerClass(CallStructuralVariationReducer.class);
		job.setNumReduceTasks(options.getReducenum());
		job.setInputFormatClass(GaeaAnySAMInputFormat.class);
		job.setMapOutputKeyClass(NewMapKey.class);
		job.setMapOutputValueClass(SamWritable.class);
		job.setOutputKeyClass(NullWritable.class);
		job.setOutputValueClass(Text.class);
		//job1.setPartitionerClass(ChrPartitionar.class);
		job.setGroupingComparatorClass(ReduceGroupingComparator.class);

		FileInputFormat.addInputPaths(job, options.getInput());
		FileOutputFormat.setOutputPath(job, new Path(options.getHdfsdir() + "/Result"));
		
		int re = job.waitForCompletion(true) ? 0:1; //submit job
		return re;
		
	}
	
	@Override
	public int run(String[] args) throws Exception {
		CallStructuralVariation sv = new CallStructuralVariation();
		int res = sv.runCallStructuralVariation(args);
		return res;
	}

}
