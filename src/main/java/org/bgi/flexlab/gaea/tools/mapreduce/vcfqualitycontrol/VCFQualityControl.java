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
package org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol;

import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.vcf.VCFMultipleInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.VCFHdfsWriter;
import org.bgi.flexlab.gaea.data.mapreduce.util.HdfsFileManager;
import org.bgi.flexlab.gaea.data.structure.header.MultipleVCFHeader;
import org.bgi.flexlab.gaea.data.structure.vcf.report.ReportDatum;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.variantrecalibratioin.VariantRecalibrationMapper;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.variantrecalibratioin.VariantRecalibrationReducer;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.HardFilter;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.VCFRecalibrator;
import org.seqdoop.hadoop_bam.KeyIgnoringVCFOutputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.IOException;

public class VCFQualityControl extends ToolsRunner{
	public static final String VQS_LOD_KEY = "VQSLOD"; // Log odds ratio of being a true variant versus being false under the trained gaussian mixture model
    public static final String CULPRIT_KEY = "culprit"; // The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out
    
    public VCFQualityControl() {
		this.toolsDescription = "Gaea VCF quality control. This module contains two options to"
				+ "perform quality control.\n"
				+ "1.Gaea variant quality score recalibration\n"
				+ "The purpose of variant recalibration is to assign a well-calibrated "
				+ "probability to each variant call in a call set.\n"
				+ "2.Gaea variant hard filter\n"
				+ "The purpose of hard filter is to filter out variants that"
				+ "don't meet the criterion of specific annotation value.\n"
				+ "Each options will automatically produce a vcf statistics report";	
	}
    
    private VCFQualityControlOptions options;
    private MultipleVCFHeader vcfHeaders;

	@Override
	public int run(String[] args) throws Exception {
		BioJob job = BioJob.getInstance();
		Configuration conf = job.getConfiguration();
		conf.setBoolean(GaeaVCFOutputFormat.HEADER_MODIFY, true);
		
		String[] remainArgs = remainArgs(args, conf);
		
		options = new VCFQualityControlOptions();
		options.parse(remainArgs);
		options.setHadoopConf(remainArgs, conf);
		
		//merge header
		vcfHeaders = new MultipleVCFHeader();
		vcfHeaders.mergeHeader(new Path(options.getInputs()), options.getOutputPath(), job, false);
				
		if(options.isRecal()) {
//			vqsr
			job.setJobName("Gaea variant quality score recalibration");
			job.setJarByClass(VCFQualityControl.class);
			job.setMapperClass(VariantRecalibrationMapper.class);
			job.setReducerClass(VariantRecalibrationReducer.class);
			job.setOutputKeyValue(IntWritable.class, Text.class, 
					NullWritable.class, VariantContextWritable.class);
			job.setNumReduceTasks(vcfHeaders.getFileNum());
			
			FileInputFormat.addInputPaths(job, options.getInputs());
			job.setInputFormatClass(VCFMultipleInputFormat.class);
	
			Path statistics = new Path(options.getOutputPath() + "/tmp");
			//job.setOutputFormatClass(VariantRecalibrationOutputFormat.class);
			job.setOutputFormatClass(GaeaVCFOutputFormat.class);
			FileOutputFormat.setOutputPath(job, statistics);
			
			return job.waitForCompletion(true) ? 0 : 1;
		} else {
//			hard filter
			HardFilter fe = new HardFilter(options.getFilterName(), options.getSnpFilter(), options.getIndelFilter());
			for(Path vcfPath : options.getInputFiles()) {
				AsciiLineReaderIterator iterator = new AsciiLineReaderIterator(new AsciiLineReader(HdfsFileManager.getInputStream(vcfPath, conf)));
				VCFCodec codec = new VCFCodec();
				VCFHeader vHeader = (VCFHeader)(((FeatureCodecHeader)codec.readHeader(iterator)).getHeaderValue());
				String output = options.getOutputPath() + vcfPath.toString();
				for(VCFHeaderLine headerLine : fe.filterHeaders())
					vHeader.addMetaDataLine(headerLine);
				
				VCFHdfsWriter writer = new VCFHdfsWriter(output, false, false, conf);
				writer.writeHeader(vHeader);
				
				ReportDatum report = null;
				while(iterator.hasNext()){
					VariantContext vc = codec.decode(iterator.next());
					vc = fe.filter(vc);
					ReportDatum datum = new ReportDatum.Builder(vc).isSnp().isIndel().isTransition().build();
					if(report == null)
						report = datum;
					else
						report.combine(datum);
					writer.add(vc);
				}
			}
			return 1;
		}
	}
}

final class VariantRecalibrationOutputFormat<K> extends FileOutputFormat<K, VariantContextWritable>{
	
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
		baseOF.setHeader(VCFRecalibrator.finalHeader);
	
		return baseOF.getRecordWriter(context, getDefaultWorkFile(context, ""));
	}

	// Allow the output directory to exist.
	@Override public void checkOutputSpecs(JobContext job) {}
}