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
package org.bgi.flexlab.gaea.tools.mapreduce.annotator;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IOUtils;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.LazyOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.mapreduce.lib.partition.InputSampler;
import org.apache.hadoop.mapreduce.lib.partition.TotalOrderPartitioner;
import org.bgi.flexlab.gaea.data.mapreduce.input.txt.AnnotationTextInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.partitioner.AnnoSortPartitioner;
import org.bgi.flexlab.gaea.data.mapreduce.partitioner.FirstPartitioner;
import org.bgi.flexlab.gaea.data.mapreduce.writable.PairWritable;
import org.bgi.flexlab.gaea.data.mapreduce.writable.VcfLineWritable;
import org.bgi.flexlab.gaea.data.structure.header.SingleVCFHeader;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.tools.annotator.config.Config;
import org.bgi.flexlab.gaea.util.FileIterator;
import org.seqdoop.hadoop_bam.VCFOutputFormat;
import org.seqdoop.hadoop_bam.util.BGZFCodec;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

public class Annotator extends ToolsRunner {

    public static final String SORT_TEMP_INFO = "sortinfo.tsv";

    private Configuration conf;
    private AnnotatorOptions options;
    private List<String> sampleNames;
    private List<String> contigs;
    private List<String> fileNames;
    private String sortInput;
    private Path sortTmpPath;

    public Annotator(){
        sampleNames = new ArrayList<>();
        fileNames = new ArrayList<>();
    }

    public Annotator(String[] arg0) throws IOException {
        sampleNames = new ArrayList<>();
        fileNames = new ArrayList<>();
        contigs = new ArrayList<>();

        conf = new Configuration();
        String[] remainArgs = remainArgs(arg0, conf);

        options = new AnnotatorOptions();
        options.parse(remainArgs);
        options.setHadoopConf(remainArgs, conf);

        sortInput = options.getInputFilePath();
        sortTmpPath = new Path(options.getTmpPath()+"/sort");
    }

    public AnnotatorOptions getOptions() {
        return options;
    }

    private int runAnnotator() throws Exception {

        Path inputPath = new Path(options.getInputFilePath());
        FileSystem fs = inputPath.getFileSystem(conf);
        FileStatus[] files = fs.listStatus(inputPath);

        for(FileStatus file : files) {//统计sample names
            System.out.println(file.getPath());
            if (file.isFile()) {
                SingleVCFHeader singleVcfHeader = new SingleVCFHeader();
                singleVcfHeader.readHeaderFrom(file.getPath(), fs);
                VCFHeader vcfHeader = singleVcfHeader.getHeader();

                fileNames.add(file.getPath().getName());

                for(String sample: vcfHeader.getSampleNamesInOrder()) {
                    if(!sampleNames.contains(sample))
                        sampleNames.add(sample);
                }
                for(SAMSequenceRecord samSequenceRecord: vcfHeader.getSequenceDictionary().getSequences()){
                    String contigName = samSequenceRecord.getSequenceName();
                    if(!contigs.contains(contigName))
                        contigs.add(contigName);
                }
            }
        }

        conf.set(VCFOutputFormat.OUTPUT_VCF_FORMAT_PROPERTY, "VCF");
        BioJob job = BioJob.getInstance(conf);

        job.setJobName("GaeaAnnotator");
        job.setJarByClass(this.getClass());
        job.setMapperClass(AnnotationMapper.class);
        job.setReducerClass(AnnotationReducer.class);
        job.setNumReduceTasks(options.getReducerNum());

        job.setMapOutputKeyClass(Text.class);
        job.setMapOutputValueClass(VcfLineWritable.class);

        job.setOutputKeyClass(Text.class);
        job.setOutputValueClass(Text.class);

        job.setInputFormatClass(TextInputFormat.class);
        FileInputFormat.addInputPath(job, inputPath);

        if(options.getRunStep() == AnnotatorOptions.RunStep.ANN){
            FileOutputFormat.setOutputPath(job, new Path(options.getOutputPath()+"/anno"));
            createSortInfo(options.getOutputPath()+"/" + SORT_TEMP_INFO);
        }else {
            Path partTmp = new Path(options.getTmpPath()+"/anno");
            FileOutputFormat.setOutputPath(job, partTmp);
            sortInput = partTmp.toString();
            createSortInfo(options.getTmpPath()+"/" + SORT_TEMP_INFO);
        }

        return job.waitForCompletion(true) ? 0 : 1;
    }

    private int runAnnoSort() throws Exception {
        if(options.getOutputFormat() == AnnotatorOptions.OutputFormat.VCF){
            return runVCFSort();
        }
        return runTSVSort();
    }

    private int runVCFSort() throws Exception {

        conf.set("io.compression.codecs", BGZFCodec.class.getCanonicalName());
        BioJob job = BioJob.getInstance(conf);

        job.setJobName("GaeaAnnotatorSort");
        job.setJarByClass(this.getClass());
        job.setMapperClass(AnnoSortMapper.class);
        job.setReducerClass(AnnoSortReducer.class);
        job.setNumReduceTasks(fileNames.size());

        job.setMapOutputKeyClass(PairWritable.class);
        job.setMapOutputValueClass(Text.class);

        job.setPartitionerClass(FirstPartitioner.class);
        job.setOutputKeyClass(NullWritable.class);
        job.setOutputValueClass(Text.class);
        job.setInputFormatClass(TextInputFormat.class);
        LazyOutputFormat.setOutputFormatClass(job, TextOutputFormat.class);
        conf.setInt("mapreduce.reduce.memory.mb", 20480);

        Path inputPath = new Path(sortInput);
        FileInputFormat.setInputPaths(job, inputPath);
        FileOutputFormat.setOutputPath(job, sortTmpPath);
        FileOutputFormat.setCompressOutput(job, true);
        FileOutputFormat.setOutputCompressorClass(job, BGZFCodec.class);

        FileSystem fs = sortTmpPath.getFileSystem(conf);
        if(job.waitForCompletion(true)){
            for (String fileName : fileNames){
                Path outputPart = getSampleOutputPath(fileName);
                Path outputName = new Path(options.getOutputPath() + "/" + fileName);
                fs.rename(outputPart, outputName);
            }
            return 0;
        }
        return 1;
    }


    private int runTSVSort() throws Exception {

        conf.set("io.compression.codecs", BGZFCodec.class.getCanonicalName());
//        conf.setInt("mapreduce.reduce.memory.mb", 30480);
        BioJob job = BioJob.getInstance(conf);

        job.setJobName("GaeaAnnotatorSort");
        job.setJarByClass(this.getClass());
        job.setMapperClass(AnnotationSortMapper.class);
        job.setReducerClass(AnnotationSortReducer.class);

        job.setMapOutputKeyClass(LongWritable.class);
        job.setOutputKeyClass(NullWritable.class);
        job.setOutputValueClass(Text.class);
        job.setInputFormatClass(AnnotationTextInputFormat.class);
        LazyOutputFormat.setOutputFormatClass(job, TextOutputFormat.class);

        if(options.getRunStep() == AnnotatorOptions.RunStep.SORT) {
            getSortInfo(job.getConfiguration(), options.getInputFilePath() + "/" + SORT_TEMP_INFO);
            FileInputFormat.setInputPaths(job, new Path(options.getInputFilePath() + "/anno"));
        }else
            FileInputFormat.setInputPaths(job, new Path(sortInput));
        FileInputFormat.setMinInputSplitSize(job, 268435456);
        FileOutputFormat.setOutputPath(job, sortTmpPath);
        FileOutputFormat.setCompressOutput(job, true);
        FileOutputFormat.setOutputCompressorClass(job, BGZFCodec.class);

        if(sampleNames.size() == 1) {
            job.setNumReduceTasks(options.getReducerNum());
            job.setPartitionerClass(TotalOrderPartitioner.class);
            Path partitionFile = new Path(options.getTmpPath() + "/_partitons.lst");
            TotalOrderPartitioner.setPartitionFile(job.getConfiguration(), partitionFile);
            System.err.println("anno-sort :: Sampling...");
            InputSampler.writePartitionFile(
                    job,
                    new InputSampler.RandomSampler<LongWritable, Text>(
                            0.01, 1000, options.getReducerNum()));
        }else {
            job.setNumReduceTasks(sampleNames.size());
            job.setPartitionerClass(AnnoSortPartitioner.class);
        }

        FileSystem fs = sortTmpPath.getFileSystem(conf);
        if(job.waitForCompletion(true)){
            Config userConfig = new Config(conf);
            Path headerPath = new Path(options.getTmpPath()+"/header.tsv.gz");
            BGZFCodec bgzfCodec = new BGZFCodec();
            OutputStream os = bgzfCodec.createOutputStream(fs.create(headerPath));
            os.write(userConfig.getHeaderString().getBytes());
            os.write('\n');
            os.close();

            for (String sample : sampleNames){
                Path outputName = new Path(options.getOutputPath() + "/" + sample + ".tsv.gz");
                OutputStream outgz = outputName.getFileSystem(conf).create(outputName);
                FSDataInputStream ins = fs.open(headerPath);
                IOUtils.copyBytes(ins, outgz, conf, false);
                FileStatus[] fileStatuses = fs.globStatus(new Path(sortTmpPath.toString() + "/" + sample + "*-r-[0-9]*"));
                for(FileStatus fstat: fileStatuses){
                    ins = fs.open(fstat.getPath());
                    IOUtils.copyBytes(ins, outgz, conf, false);
                    ins.close();
                }
                outgz.write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
                outgz.close();
            }
            fs.delete(new Path(options.getTmpPath()), true);
            return 0;
        }
        return 1;
    }

    private Path getSampleOutputPath(String sample) throws IOException {
        Path outputPath = new Path(options.getOutputPath());
        FileSystem fs = outputPath.getFileSystem(conf);
        FileStatus[] fileStatuses = fs.globStatus(new Path(sortTmpPath.toString() + "/" + sample + "-r-[0-9]*"));
        if(fileStatuses.length == 0){
            System.err.println(sample+": cann't get the output part file!");
            FileStatus[] fss = fs.globStatus(new Path(options.getOutputPath() + "/*"));
            for (FileStatus f: fss){
                System.err.println("OutPath" + f.getPath().toString());
            }
            return null;
        }
        return fileStatuses[0].getPath();
    }

    private void createSortInfo(String sortInfo) throws IOException {
        conf.setStrings("sampleName", sampleNames.toArray(new String[0]));
        conf.setStrings("contigName", contigs.toArray(new String[0]));

        Path p = new Path(sortInfo);
        OutputStream out = p.getFileSystem(conf).create(p);
        out.write(String.join("\t", sampleNames).getBytes());
        out.write('\n');
        out.write(String.join("\t", contigs).getBytes());
        out.write('\n');
        out.close();
    }

    private void getSortInfo(Configuration conf, String sortInfo) throws IOException {
        boolean firstLine = true;
        FileIterator fi = new FileIterator(sortInfo);
        while (fi.hasNext()){
            String line = fi.next().toString();
            if(firstLine) {
                String[] samples = line.split("\t");
                conf.setStrings("sampleName", samples);
                for (String sample: samples){
                    if(!sampleNames.contains(sample))
                        sampleNames.add(sample);
                }
                firstLine = false;
            }else {
                String[] fields = line.split("\t");
                conf.setStrings("contigName", fields);
            }
        }
        fi.close();
    }

    @Override
    public int run(String[] args) throws Exception {
        Annotator annotator = new Annotator(args);
        if(annotator.getOptions().getRunStep() == AnnotatorOptions.RunStep.SORT)
            return annotator.runAnnoSort();
        else if(annotator.getOptions().getRunStep() == AnnotatorOptions.RunStep.ANN)
            return annotator.runAnnotator() == 0 ? 0 : 1;

        return annotator.runAnnotator() == 0 ? annotator.runAnnoSort() : 1;
    }

}
