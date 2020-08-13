package org.bgi.flexlab.gaea.tools.mapreduce.combinegvcfs;

import java.util.Set;
import java.util.TreeSet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.VCFHdfsWriter;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCalling;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCallingMapper;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCallingOptions;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCallingReducer;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCallingVCFInputFormat;
import org.bgi.flexlab.gaea.util.Utils;
import org.seqdoop.hadoop_bam.KeyIgnoringVCFOutputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.*;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;

public class CombineGVCFs extends ToolsRunner{
	public final static String INPUT_ORDER = "input.name.order";
	public CombineGVCFs(){
		this.toolsDescription = "Combine for gvcfs";
	}
	private Set<String> getSampleList(Set<VCFHeader> headers){
		Set<String> samples = new TreeSet<String>();
		for(VCFHeader header : headers){
			for ( String sample : header.getGenotypeSamples() ) {
				samples.add(GaeaGvcfVariantContextUtils.mergedSampleName(null, sample, false));
			}
		}
		
		return samples;
	}
	private VCFHeader getVCFHeaderFromInput(Set<VCFHeader> headers) {
		Set<String> samples = getSampleList(headers);
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(headers, true);
        VCFHeader vcfHeader = new VCFHeader(headerLines, samples);
        
        headers.clear();
        samples.clear();
        headerLines.clear();
        
        return vcfHeader;
	}
	@Override
	public int run(String[] args) throws Exception {
		BioJob job = BioJob.getInstance();
        Configuration conf = job.getConfiguration();
        
        String[] remainArgs = remainArgs(args, conf);
        CombineOptions options = new CombineOptions();
        options.parse(remainArgs);
        options.setHadoopConf(remainArgs, conf);
        conf.set(KeyIgnoringVCFOutputFormat.OUTPUT_VCF_FORMAT_PROPERTY, options.getOuptputFormat().toString());
        conf.setBoolean(GaeaVCFOutputFormat.HEADER_MODIFY, true);
        MultipleVCFHeaderForJointCalling multiVcfHeader = new MultipleVCFHeaderForJointCalling();
        if(options.getVcfHeaderFile() != null) {
            conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVcfHeaderFile());
//            multiVcfHeader.headersConfig(new Path(options.getVcfHeaderFile()), options.getVCFHeaderOutput()+"/vcfHeaders", conf);
        }else {
            conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVCFHeaderOutput() + "/vcfFileHeader.vcf");
            multiVcfHeader.headersConfig(options.getInput(), options.getVCFHeaderOutput()+"/vcfHeaders", conf);
//            VCFHeader vcfHeader = getVCFHeaderFromInput(multiVcfHeader.getHeaders());
            VCFHdfsWriter vcfHdfsWriter = new VCFHdfsWriter(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP), false, false, conf);
//            vcfHdfsWriter.writeHeader(vcfHeader);
            vcfHdfsWriter.close();
        }
//        conf.set(INPUT_ORDER, Utils.join(",", multiVcfHeader.getSamplesAsInputOrder()));

        job.setJobName("Gaea CombineGVCFs");
        
        job.setJarByClass(CombineGVCFs.class);
        job.setWindowsBasicMapperClass(JointCallingMapper.class, options.getWindowsSize(),0);
        job.setReducerClass(CombineGVCFsReducer.class);
        //job.setCombinerClass(JointCallingPrepareReducer.class);
        
        job.setNumReduceTasks(options.getReducerNumber());
        job.setOutputKeyValue(WindowsBasedWritable.class,VariantContextWritable.class, NullWritable.class, VariantContextWritable.class);
        
        job.setInputFormatClass(JointCallingVCFInputFormat.class);
		job.setOutputFormatClass(GaeaVCFOutputFormat.class);
        
        FileInputFormat.setInputPaths(job, options.getInput().toArray(new Path[options.getInput().size()]));
		FileOutputFormat.setOutputPath(job, new Path(options.getOutput()));
		
		return job.waitForCompletion(true) ? 0 : 1;
	}
	
	
}
