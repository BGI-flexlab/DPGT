package org.bgi.flexlab.gaea.data.datawarehouse.genomics;

import java.io.IOException;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.HColumnDescriptor;
import org.apache.hadoop.hbase.HTableDescriptor;
import org.apache.hadoop.hbase.client.HBaseAdmin;
import org.apache.hadoop.hbase.client.HTable;
import org.apache.hadoop.hbase.client.Put;
import org.apache.hadoop.hbase.io.ImmutableBytesWritable;
import org.apache.hadoop.hbase.mapreduce.HFileOutputFormat2;
import org.apache.hadoop.hbase.mapreduce.LoadIncrementalHFiles;
import org.apache.hadoop.hbase.mapreduce.PutSortReducer;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.vcf.VCFMultipleInputFormat;
import org.bgi.flexlab.gaea.data.structure.header.MultipleVCFHeader;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;

public class LoadVCFToHBase extends ToolsRunner{
	
	private LoadVCFToHBaseOptions options = null;
	
	private void createTable(Configuration conf, String tableName) throws IOException {
		HBaseAdmin admin = new HBaseAdmin(conf);
		if (admin.tableExists(tableName)) {
			admin.close();
			return;
		}
		HTableDescriptor table = new HTableDescriptor(tableName);
		table.addFamily(new HColumnDescriptor("info"));
		admin.createTable(table);
		admin.close();
	}
	
	private void LoadHFile2HBase(Configuration conf,String tableName,String hfile) throws Exception{
		conf.set("hbase.metrics.showTableName", "false");
		LoadIncrementalHFiles loader = new LoadIncrementalHFiles(conf);
		HBaseAdmin admin = new HBaseAdmin(conf);
		HTable table = new HTable(conf, tableName);

		loader.doBulkLoad(new Path(hfile), table);
		table.flushCommits();
		table.close();
		admin.close();
	}

	@Override
	public int run(String[] args) throws Exception {
		Configuration conf = HBaseConfiguration.create();
		
		String[] remainArgs = remainArgs(args, conf);
		options = new LoadVCFToHBaseOptions();
		options.parse(remainArgs);
		options.setHadoopConf(remainArgs, conf);
		
		conf.addResource(new Path(options.getConfig() + "hbase-site.xml"));
    	conf.addResource(new Path(options.getConfig() + "core-site.xml"));
    	conf.set("vcfHeader", options.getHeaderOutput());
		Job job = new Job(conf);

		createTable(conf,options.getTableName());
		
		MultipleVCFHeader vcfHeaders = new MultipleVCFHeader();
		vcfHeaders.mergeHeader(new Path(options.getInput()),options.getHeaderOutput(), job, false);

		job.setJobName("vcf to hbase");
		job.setNumReduceTasks(options.getReducerNumber());
		job.setInputFormatClass(VCFMultipleInputFormat.class);

		job.setJarByClass(LoadVCFToHBase.class);
		job.setMapperClass(VCFToHBaseMapper.class);
		job.setReducerClass(PutSortReducer.class);

		job.setMapOutputKeyClass(ImmutableBytesWritable.class);
		job.setMapOutputValueClass(Put.class);
		job.setOutputKeyClass(ImmutableBytesWritable.class);
		job.setOutputValueClass(Put.class);

		FileInputFormat.setInputPaths(job, new Path(options.getInput()));
		FileOutputFormat.setOutputPath(job, new Path(options.getHFileOutput()));
		
		HFileOutputFormat2.configureIncrementalLoad(job, new HTable(conf,options.getTableName()));

		if (job.waitForCompletion(true)) {
			LoadHFile2HBase(conf,options.getTableName(),options.getHFileOutput());
			return 0;
		} else {
			return 1;
		}
	}

}
