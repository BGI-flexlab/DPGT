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
package org.bgi.flexlab.gaea.tools.mapreduce.bamsort;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.Log;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IOUtils;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.LazyOutputFormat;
import org.apache.hadoop.mapreduce.lib.partition.InputSampler;
import org.apache.hadoop.mapreduce.lib.partition.TotalOrderPartitioner;
import org.bgi.flexlab.gaea.data.mapreduce.input.bam.GaeaAnySAMSortInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.bam.GaeaNoheaderBamOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.mapreduce.input.bam.GaeaSamSortRecordReader;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.seqdoop.hadoop_bam.SAMFormat;
import org.seqdoop.hadoop_bam.util.SAMOutputPreparer;
import org.seqdoop.hadoop_bam.util.Timer;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BamSort extends ToolsRunner {

    BamSortOptions options;
    final String intermediateOutName = "hadoopBamSort";
    private SAMFileHeader header;
    private Map<String, String> formatSampleName;
    private List<String> sampleNames;
    private List<Path> inputs = new ArrayList<>();
    BioJob job;
    Configuration conf;
    Path tmpPath;
    Path tmpDir;
    Path outputPath;
    SAMFormat format;
    String fileSuffix;

    public BamSort(){
        this.toolsDescription = "Gaea BamSort";
        sampleNames = new ArrayList<>();
    }


    public int runSort(String[] args) throws IOException, ClassNotFoundException, InterruptedException, URISyntaxException {
        job = BioJob.getInstance();
        conf = job.getConfiguration();
        String[] remainArgs = remainArgs(args, conf);
        options = new BamSortOptions();
        options.parse(remainArgs);
        tmpPath = new Path(options.getTmpPath()+"/data");
        tmpDir = new Path(options.getTmpPath());
        outputPath = new Path(options.getOutdir());
        options.setHadoopConf(remainArgs, conf);

        inputs = options.getInputs();
        header = job.setHeader(inputs, outputPath);

        format = options.getOutputFormat(conf);
        fileSuffix = BamSortUtils.getFileSuffix(format);

        for (SAMReadGroupRecord rg : header.getReadGroups()) {
            if(!sampleNames.contains(rg.getSample()))
                sampleNames.add(rg.getSample());
        }
        conf.setStrings(GaeaSamSortRecordReader.SAMPLENAME_ARRAY_PROP, sampleNames.toArray(new String[0]));

        conf.setBoolean(SortOutputFormat.WRITE_HEADER_PROP, false);
        conf.set(SortOutputFormat.OUTPUT_NAME_PROP, intermediateOutName);
        conf.set(SortOutputFormat.OUTPUT_SAM_FORMAT_PROPERTY, options.getOutputFormat());


        job.setJarByClass(BamSort.class);
        job.setMapperClass(Mapper.class);
        job.setReducerClass(BamSortReducer.class);
        job.setJobName("bamsort");
        job.setNumReduceTasks(options.getReducerNum());
        job.setMapOutputKeyClass(LongWritable.class);
        job.setOutputKeyClass(NullWritable.class);
        job.setOutputValueClass(SamRecordWritable.class);
        job.setInputFormatClass(GaeaAnySAMSortInputFormat.class);
        if (format == SAMFormat.BAM)
            LazyOutputFormat.setOutputFormatClass(job, GaeaNoheaderBamOutputFormat.class);

        for (Path in : inputs)
            FileInputFormat.addInputPath(job, in);

        FileOutputFormat.setOutputPath(job, tmpPath);

        job.setPartitionerClass(TotalOrderPartitioner.class);

        Path partitionFile;
        if(options.getPartitionFile() == null) {
            partitionFile = new Path(tmpDir + "/_partitons.lst");
            TotalOrderPartitioner.setPartitionFile(job.getConfiguration(), partitionFile);
            System.err.println("bam-sort :: Sampling...");
            InputSampler.writePartitionFile(
                    job,
                    new InputSampler.RandomSampler<LongWritable, SamRecordWritable>(
                            0.01, 10000, options.getReducerNum()));

        }else {
            System.err.println("bamsort :: use partitionFile:"+options.getPartitionFile() + " ...");
            partitionFile = new Path(options.getPartitionFile());
            TotalOrderPartitioner.setPartitionFile(job.getConfiguration(), partitionFile);
        }


        if(job.waitForCompletion(true)) {
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

            FileSystem fs = tmpPath.getFileSystem(conf);
            for (String sample : sampleNames){
                SAMFileHeader sampleHeader = BamSortUtils.deleteSampleFromHeader(header, sample);
                Path outputPath = new Path(options.getOutdir() + "/" + sample + fileSuffix);
                final FileSystem dstFS = outputPath.getFileSystem(conf);

                OutputStream os = dstFS.create(outputPath);
                new SAMOutputPreparer().prepareForRecords(
                        os, format,
                        sampleHeader);

                FileStatus[] fileStatuses = fs.globStatus(new Path(tmpPath.toString() + "/" + sample + "*-r-[0-9]*"));
                for(FileStatus fstat: fileStatuses){
                    FSDataInputStream ins = fs.open(fstat.getPath());
                    IOUtils.copyBytes(ins, os, conf, false);
                    ins.close();
                }
                os.write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
                os.close();
            }
//            fs.delete(tmpPath, true);
            return 0;
        }
        return 1;
    }

    public int runMultiSort2() throws IOException, ClassNotFoundException, InterruptedException, URISyntaxException {

        conf.setBoolean(SortOutputFormat.WRITE_HEADER_PROP, options.getOutdir() == null);
        conf.set(SortOutputFormat.OUTPUT_NAME_PROP, intermediateOutName);
        conf.set(SortOutputFormat.OUTPUT_SAM_FORMAT_PROPERTY, options.getOutputFormat());

        job.setMapperClass(SortMapper.class);
        job.setReducerClass(SortReducer.class);
        job.setJobName("multi bamsort2");

        job.setMapOutputKeyClass(LongWritable.class);
        job.setOutputKeyClass(NullWritable.class);
        job.setOutputValueClass(SamRecordWritable.class);

        job.setInputFormatClass(GaeaAnySAMSortInputFormat.class);
        job.setOutputFormatClass(SortOutputFormat.class);

        Timer t = new Timer();
        formatSampleName = new HashMap<>();

        for (final Path in : inputs)
            FileInputFormat.addInputPath(job, in);

        FileOutputFormat.setOutputPath(job, tmpPath);
        if(options.getRenames() != null){
            header = BamSortUtils.replaceSampleName(header.clone(), options.getRenames());
        }

        for (SAMReadGroupRecord rg : header.getReadGroups()) {
            String fsn = BamSortUtils.formatSampleName(rg.getSample());
            if (!formatSampleName.containsKey(fsn)) {
                formatSampleName.put(fsn, rg.getSample());
                SortMultiOutputs.addNamedOutput(job, fsn,
                        SortOutputFormat.class, NullWritable.class,
                        SamRecordWritable.class);
            }
        }

        job.setPartitionerClass(SamSortPartition.class);

        t.start();

        if (!job.waitForCompletion(true)) {
            System.err.println("sort :: Job failed.");
            return 4;
        }

        System.out.printf("sort :: Job complete in %d.%03d s.\n",
                t.stopS(), t.fms());

        return 0;
    }

    public int mergeBAM(Configuration conf, SAMFormat format) {
        Timer t = new Timer();
        try {
            Log.setGlobalLogLevel(Log.LogLevel.ERROR);
            System.out.println("sort :: Merging output...");


            final FileSystem srcFS = tmpPath.getFileSystem(conf);
            FileSystem dstFS = outputPath.getFileSystem(conf);

            // create output stream foreach sample
            Map<String, OutputStream> outs = new HashMap<>();
            String fileSuffix = BamSortUtils.getFileSuffix(format);
            for (String sampleName : sampleNames) {
                Path sPath = new Path(options.getOutdir() + "/" + sampleName + fileSuffix);
                OutputStream os = dstFS.create(sPath);
                outs.put(sampleName, os);
            }

            // Then, the actual SAM or BAM contents.
            for (String fsn : formatSampleName.keySet()) {
                t.start();
                SAMFileHeader newHeader = BamSortUtils.deleteSampleFromHeader(header,formatSampleName
                        .get(fsn));
                new SAMOutputPreparer().prepareForRecords(
                        outs.get(formatSampleName.get(fsn)), format,
                        newHeader);
                final FileStatus[] parts = srcFS.globStatus(new Path(
                        options.getTmpPath(), fsn + "-*-[0-9][0-9][0-9][0-9][0-9]*"));

                for (final FileStatus part : parts) {
                    System.out.flush();

                    final InputStream ins = srcFS.open(part.getPath());
                    IOUtils.copyBytes(ins,
                            outs.get(formatSampleName.get(fsn)), conf,
                            false);
                    ins.close();
                }

                if (format == SAMFormat.BAM)
                    outs.get(formatSampleName.get(fsn))
                            .write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
                outs.get(formatSampleName.get(fsn)).close();
                System.out.printf(
                        "sort :: Merging " + formatSampleName.get(fsn)
                                + " complete in %d.%03d s.\n", t.stopS(),
                        t.fms());
            }
        } catch (IOException e) {
            System.err.printf("sort :: Output merging failed: %s\n", e);
            return 5;
        }
        return 0;
    }


    private Path getSampleOutputPath(String sample) throws IOException {
        FileSystem fs = tmpPath.getFileSystem(conf);
        FileStatus[] fileStatuses = fs.globStatus(new Path(options.getTmpPath() + "/" + sample + "-r-[0-9]*"));
        if(fileStatuses.length == 0){
            System.err.println(sample+": cann't get the output part file!");
            FileStatus[] fss = fs.globStatus(new Path(options.getTmpPath() + "/*"));
            for (FileStatus f: fss){
                System.err.println("OutPath" + f.getPath().toString());
            }
            return null;
        }
        return fileStatuses[0].getPath();
    }

    public static void configureSampling(
            Path workDir, String outName, Configuration conf)
            throws IOException
    {
        final Path partition =
                workDir.getFileSystem(conf).makeQualified(
                        new Path(workDir, "_partitioning" + outName));

        TotalOrderPartitioner.setPartitionFile(conf, partition);
        try {
            final URI partitionURI = new URI(
                    partition.toString() + "#" + partition.getName());

            if (partitionURI.getScheme().equals("file"))
                return;

            DistributedCache.addCacheFile(partitionURI, conf);
            DistributedCache.createSymlink(conf);
        } catch (URISyntaxException e) { throw new RuntimeException(e); }
    }

    @Override
    public int run(String[] args) throws Exception {
        BamSort sort = new BamSort();
        return sort.runSort(args);
    }
}