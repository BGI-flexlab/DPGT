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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2010 Aalto University 
 *
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.mapreduce.vcf.sort;


import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.*;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.MultipleOutputs;
import org.apache.hadoop.mapreduce.lib.partition.InputSampler;
import org.apache.hadoop.mapreduce.lib.partition.TotalOrderPartitioner;
import org.bgi.flexlab.gaea.data.mapreduce.input.vcf.VCFRecordReader;
import org.bgi.flexlab.gaea.data.mapreduce.util.HdfsFileManager;
import org.bgi.flexlab.gaea.data.structure.header.MultipleVCFHeader;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.util.SortUilts;
import org.seqdoop.hadoop_bam.KeyIgnoringVCFOutputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.HashMap;
import java.util.Map;

public class VCFSort extends ToolsRunner{

	private static VCFSortOptions options;
	
	public VCFSort() {
		this.toolsDescription = "Gaea vcf sorting\n"
				+ "a vcf sorting tool that can process multiple vcf files a time";
	}
	
	@Override
	public int run(String[] args) {
		try {
			BioJob job = BioJob.getInstance();
			Configuration conf = job.getConfiguration();
			String[] remainArgs = remainArgs(args, conf);
			options = new VCFSortOptions();
			options.parse(remainArgs);

			conf.set(SortOutputFormat.INPUT_PATH_PROP, options.getInputFileList().get(0).toString());

			SortUilts.configureSampling(new Path(options.getTempOutput()), job, options);

			MultipleVCFHeader mVcfHeader = mergeHeader(options, job);

			initChrOrder(conf);

			options.setHadoopConf(remainArgs, conf);

			job.setJobName("VCFSort");
			
			job.setJarByClass(VCFSort.class);
			job.setMapperClass(Mapper.class);;
			job.setReducerClass(VCFSortReducer.class);
			
			job.setNumReduceTasks(options.getReducerNumber());
	
			job.setPartitionerClass(TotalOrderPartitioner.class);
			
			job.setOutputKeyValue(LongWritable.class, VariantContextWritable.class,
					NullWritable.class,VariantContextWritable.class);
			job.setInputFormatClass(SortInputFormat.class);
			job.setOutputFormatClass(SortOutputFormat.class);
			
			FileInputFormat.setInputPaths(job, 
					options.getInputFileList().toArray(new Path[options.getInputFileList().size()]));
			FileOutputFormat.setOutputPath(job, new Path(options.getWorkPath()));
	
			setMultiOutputs(mVcfHeader, job);
			
			InputSampler.<LongWritable, VariantContextWritable>writePartitionFile(
					job,
					new InputSampler.RandomSampler<LongWritable,VariantContextWritable>
						(0.01, 300, options.getReducerNumber()));
			
			job.submit();
			
			if (!job.waitForCompletion(true)) {
				System.err.println("vcf-Sort :: Job failed.");
				return 4;
			}
			
			if (options.getOutputPath() != null) {
				SortUilts.merge(mVcfHeader,options, conf);
			}
		} catch (IOException e) {
			System.err.printf("vcf-Sort :: Hadoop error: %s\n", e);
			return 4;
		} catch (ClassNotFoundException e) { 
			throw new RuntimeException(e); 
		} catch (InterruptedException e) { 
			throw new RuntimeException(e); 
		}
	
		return 0;
	}
	
	private MultipleVCFHeader mergeHeader(VCFSortOptions options, Job job) {
		MultipleVCFHeader mVcfHeader = new MultipleVCFHeader();
		mVcfHeader.mergeHeader(new Path(options.getInput()), options.getOutputPath(), job, false);
		return mVcfHeader;
	}
	
	private void initChrOrder(Configuration conf) throws IOException {
		Map<String, Long> chrOrder = new HashMap<>();
		FSDataInputStream ins = HdfsFileManager.getInputStream(new Path(options.getChrOrderFile()), conf);
        AsciiLineReaderIterator it = new AsciiLineReaderIterator(new AsciiLineReader(ins));
        String line;  
        long i = 1;
        while(it.hasNext()){
        	line = it.next();
        	//System.err.println("fai:" + line);
        	String[] cols = line.split("\t");
        	chrOrder.put(cols[0].trim(), i++ << 32 );
        }
        it.close();
        storage(chrOrder, conf);
	}
	
	private void storage(Map<String, Long> chrOrder, Configuration conf) throws IOException {
		ObjectOutputStream oStream = new ObjectOutputStream(HdfsFileManager.getOutputStream(
				new Path(conf.get(VCFRecordReader.CHR_ORDER_PROPERTY)), conf));
		oStream.writeObject(chrOrder);
		oStream.close();
	}
	
	private void setMultiOutputs(MultipleVCFHeader mVcfHeader, BioJob job) {
		// TODO Auto-generated method stub
		int i = 0;
		Map<Integer, String> multiOutputs = new HashMap<>();
		for(int id : mVcfHeader.getFileName2ID().values()) {
			multiOutputs.put(id, "SortResult" + ++i);
			MultipleOutputs.addNamedOutput(job, multiOutputs.get(id), SortOutputFormat.class, NullWritable.class, VariantContextWritable.class);
		}
		options.setMultiOutputs(multiOutputs);
	}
	
}

final class VCFSortReducer extends Reducer<LongWritable, VariantContextWritable, NullWritable, VariantContextWritable> {
	private MultipleOutputs<NullWritable, VariantContextWritable> mos;
	private VCFSortOptions options;
	private Map<Integer, String> multiOutputs;

	@Override 
	protected void setup(Context context) throws IOException{
		Configuration conf = context.getConfiguration();
		options = new VCFSortOptions();
		options.getOptionsFromHadoopConf(conf);
		multiOutputs = options.getMultiOutputs();
		mos = new MultipleOutputs<NullWritable, VariantContextWritable>(context);	
	}

	@Override 
	protected void reduce(
			LongWritable key, Iterable<VariantContextWritable> records,
			Reducer<LongWritable,VariantContextWritable,
			        NullWritable,VariantContextWritable>.Context
				ctx)
		throws IOException, InterruptedException
	{
		int id = (int)(key.get() >> 40);
		for (VariantContextWritable rec : records)
			//NullWirtable.get()输出null;
			mos.write(multiOutputs.get(id), NullWritable.get(), rec);
	}

	@Override
	protected void cleanup(Context context) throws IOException, InterruptedException {
		mos.close();
	}

}

class SortInputFormat extends FileInputFormat<LongWritable, VariantContextWritable> {
	
	@Override
	public RecordReader<LongWritable, VariantContextWritable> createRecordReader(InputSplit split, TaskAttemptContext ctx) throws IOException, InterruptedException {
		Configuration conf = ctx.getConfiguration();
		return new VCFRecordReader(conf, conf.get(VCFRecordReader.CHR_ORDER_PROPERTY), true);
	}
}

final class SortOutputFormat<K> extends FileOutputFormat<K, VariantContextWritable> {
	public static final String INPUT_PATH_PROP = "hadoopbam.vcfsort.inpath";
	
	private KeyIgnoringVCFOutputFormat<K> baseOF;
	
	private void initBaseOF(Configuration conf) {
		if (baseOF == null)
			baseOF = new KeyIgnoringVCFOutputFormat<K>(conf);
	}

	@Override public RecordWriter<K,VariantContextWritable> getRecordWriter(
			TaskAttemptContext context)
		throws IOException {
		final Configuration conf = context.getConfiguration();
		initBaseOF(conf);
		if (baseOF.getHeader() == null) {
			final Path p = new Path(conf.get(INPUT_PATH_PROP));
			baseOF.readHeaderFrom(p, p.getFileSystem(conf));
		}
	
		return baseOF.getRecordWriter(context, getDefaultWorkFile(context, ""));
	}

	// Allow the output directory to exist.
	@Override public void checkOutputSpecs(JobContext job) {}
}
