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
package org.bgi.flexlab.gaea.tools.mapreduce.vcfstats;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.tools.vcfstats.report.VCFReport;
import org.seqdoop.hadoop_bam.VCFInputFormat;
import org.seqdoop.hadoop_bam.VCFOutputFormat;

import java.util.concurrent.TimeUnit;

public class VCFStats extends ToolsRunner {

    private Configuration conf;
    private VCFStatsOptions options;

    private int runVCFStats(String[] arg0) throws Exception {

        conf = new Configuration();
        String[] remainArgs = remainArgs(arg0, conf);

        options = new VCFStatsOptions();
        options.parse(remainArgs);
        options.setHadoopConf(remainArgs, conf);
        conf.set(VCFOutputFormat.OUTPUT_VCF_FORMAT_PROPERTY, "VCF");
        BioJob job = BioJob.getInstance(conf);

        job.setJobName("GaeaVCFStats");
        job.setJarByClass(this.getClass());
        job.setMapperClass(VCFStatsMapper.class);
//        job.setReducerClass(VCFStatsReducer.class);
        job.setNumReduceTasks(0);

//        job.setMapOutputKeyClass(Text.class);
//        job.setMapOutputValueClass(VariantContextWritable.class);


        VCFInputFormat.addInputPath(job, new Path(options.getInput()));

        job.setOutputKeyClass(NullWritable.class);
        job.setOutputValueClass(Text.class);
//        job.setOutputKeyValue();

        job.setInputFormatClass(VCFInputFormat.class);
        job.setOutputFormatClass(TextOutputFormat.class);

        Path partTmp = new Path(options.getOutputPath() + "/out_vcf");
        FileOutputFormat.setOutputPath(job, partTmp);

//        MultipleOutputs.addNamedOutput(job, "Statistic",
//                TextOutputFormat.class, NullWritable.class, Text.class);

        if (job.waitForCompletion(true)) {
            int loop = 0;
            while (!partTmp.getFileSystem(conf).exists(partTmp) && loop < 30){
                TimeUnit.MILLISECONDS.sleep(1000);
                loop ++;
            }

            VCFReport report = new VCFReport(options);
            report.mergeReport(partTmp, conf,
                    new Path(options.getOutputPath()));
            return 0;
        } else {
            return 1;
        }
    }

    @Override
    public int run(String[] args) throws Exception {
        VCFStats vcfStats = new VCFStats();
        return vcfStats.runVCFStats(args);
    }

}
