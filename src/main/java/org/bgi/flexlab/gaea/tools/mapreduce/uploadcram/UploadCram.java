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
package org.bgi.flexlab.gaea.tools.mapreduce.uploadcram;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.cram.GaeaCramInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.bam.GaeaBamOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;

import static org.bgi.flexlab.gaea.data.mapreduce.input.cram.GaeaCramRecordReader.INPUTFORMAT_REFERENCE;

public class UploadCram  extends ToolsRunner {

    public UploadCram() {
        this.toolsDescription = "Upload cram to hdfs(in bam format)";
    }

    private UploadCramOptions options;

    @Override
    public int run(String[] args) throws Exception {
        BioJob job = BioJob.getInstance();
        Configuration conf = job.getConfiguration();
        String[] remainArgs =  remainArgs(args,conf);
        options = new UploadCramOptions();
        options.parse(remainArgs);
        options.setHadoopConf(remainArgs, conf);
        conf.set(INPUTFORMAT_REFERENCE, options.getLocalReferenceSequencePath());

        job.setJobName("UploadCram");
        job.setJarByClass(UploadCram.class);
        job.setMapperClass(UploadCramMapper.class);
        job.setOutputKeyClass(NullWritable.class);
        job.setOutputValueClass(SamRecordWritable.class);
        job.setNumReduceTasks(0);

        FileInputFormat.setInputPaths(job, options.getInputs().toArray(new Path[options.getInputs().size()]));
        job.setInputFormatClass(GaeaCramInputFormat.class);
        job.setOutputFormatClass(GaeaBamOutputFormat.class);

        FileOutputFormat.setOutputPath(job, new Path(options.getOutputPath()));
        return job.waitForCompletion(true) ? 0 : 1;
    }

    public static void main(String[] args) throws Exception {
        UploadCram uploadCram = new UploadCram();
        uploadCram.run(args);
    }
}
