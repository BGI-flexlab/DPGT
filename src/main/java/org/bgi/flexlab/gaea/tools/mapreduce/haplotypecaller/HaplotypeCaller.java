package org.bgi.flexlab.gaea.tools.mapreduce.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.VCFHdfsWriter;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.data.structure.bam.filter.HaplotypeCallerFilter;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.WindowsBasedPlaceholderSamRecordMapper;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.WindowsBasedSamRecordMapper;
import org.bgi.flexlab.gaea.tools.haplotypecaller.HaplotypeCallerTraversal;
import org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection.HaplotypeCallerArgumentCollection;
import org.bgi.flexlab.gaea.tools.haplotypecaller.writer.GVCFHadoopWriter;
import org.bgi.flexlab.gaea.util.ReadUtils;
import org.seqdoop.hadoop_bam.VCFOutputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.util.Set;

import static org.bgi.flexlab.gaea.framework.tools.mapreduce.WindowsBasedMapper.BASERECALIBRATOR_ONLY;
import static org.bgi.flexlab.gaea.framework.tools.mapreduce.WindowsBasedPlaceholderSamRecordMapper.WINDOWS_OUTPUT_ALL;

public class HaplotypeCaller extends ToolsRunner {

	@Override
	public int run(String[] args) throws Exception {
		BioJob job = BioJob.getInstance();
        Configuration conf = job.getConfiguration();
        
        String[] remainArgs = remainArgs(args, conf);
        HaplotypeCallerOptions options = new HaplotypeCallerOptions();
        options.parse(remainArgs);
        options.setHadoopConf(remainArgs, conf);
        conf.set(VCFOutputFormat.OUTPUT_VCF_FORMAT_PROPERTY, "VCF");
        conf.setBoolean(GaeaVCFOutputFormat.HEADER_MODIFY, true);
        conf.getBoolean(BASERECALIBRATOR_ONLY,true);
        conf.getBoolean(WINDOWS_OUTPUT_ALL,options.isOutputAllWindows());  //if true, output N or uncovor region windows


        SAMFileHeader samFileHeader = job.setHeader(options.getInput(), new Path(options.getHeaderOutput()));
        // /vcf header
        conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getHeaderOutput() + "/vcfFileHeader.vcf");
        HaplotypeCallerTraversal haplotypecaller = new HaplotypeCallerTraversal(null, options, samFileHeader);
        VCFHdfsWriter vcfHdfsWriter = new VCFHdfsWriter(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP), false, false, conf);
        VCFHeader vcfHeader = haplotypecaller.getVCFHeader();

        HaplotypeCallerArgumentCollection hcArgs = options.getHaplotypeCallerArguments();
        if(options.isGVCF())
            vcfHeader = new GVCFHadoopWriter(null, vcfHeader, hcArgs.GVCFGQBands, hcArgs.samplePloidy).getGVCFHeader();

        vcfHdfsWriter.writeHeader(vcfHeader);
        vcfHdfsWriter.close();

        job.setJobName("Gaea haplotype caller");
        job.setFilterClass(HaplotypeCallerFilter.class);
        
        job.setJarByClass(HaplotypeCaller.class);
        if(options.isGVCF()) {
            Set<String> sampleSet = ReadUtils.getSamplesFromHeader(samFileHeader);
//            job.setWindowsBasicMapperClass(WindowsBasedSamRecordMapper.class, options.getWindowSize(), options.getWindowsExtendSize(), sampleSet.size() > 1);
            job.setWindowsBasicMapperClass(WindowsBasedPlaceholderSamRecordMapper.class, options.getWindowSize(), options.getWindowsExtendSize(), sampleSet.size() > 1);
        }else
            job.setWindowsBasicMapperClass(WindowsBasedSamRecordMapper.class, options.getWindowSize(), options.getWindowsExtendSize());
        job.setReducerClass(HaplotypeCallerReducer.class);
        
        job.setNumReduceTasks(options.getReducerNumber());
        job.setOutputKeyValue(WindowsBasedWritable.class,SamRecordWritable.class, NullWritable.class, VariantContextWritable.class);
        
        job.setAnySamInputFormat(options.getInputFormat());
		job.setOutputFormatClass(GaeaVCFOutputFormat.class);
        
        FileInputFormat.setInputPaths(job, options.getInput().toArray(new Path[options.getInput().size()]));
		FileOutputFormat.setOutputPath(job, new Path(options.getVCFOutput()));
		
		return job.waitForCompletion(true) ? 0 : 1;
	}

}
