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
package org.bgi.flexlab.gaea.tools.mapreduce.genotyper;

import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.VCFHdfsWriter;
import org.bgi.flexlab.gaea.data.mapreduce.writable.AlignmentBasicWritable;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.data.structure.bam.filter.GenotyperFilter;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.WindowsBasedAlignmentMapper;
import org.bgi.flexlab.gaea.tools.genotyer.VariantCallingEngine;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.VariantAnnotatorEngine;
import org.seqdoop.hadoop_bam.KeyIgnoringVCFOutputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.IOException;

import static org.bgi.flexlab.gaea.framework.tools.mapreduce.WindowsBasedMapper.REFERENCE_REGION;

/**
 * Created by zhangyong on 2017/3/3.
 */
public class Genotyper extends ToolsRunner {

    public Genotyper() {
        this.toolsDescription = "Gaea genotyper";
    }

    private GenotyperOptions options = null;

    private int runGenoytper(String[] args) throws IOException, ClassNotFoundException, InterruptedException {
        BioJob job = BioJob.getInstance();
        Configuration conf = job.getConfiguration();
        String[] remainArgs = remainArgs(args, conf);

        options = new GenotyperOptions();
        options.parse(remainArgs);
        options.setHadoopConf(remainArgs, conf);
        // merge header and set to configuration
        job.setHeader(new Path(options.getInput()), new Path(options.getBAMHeaderOutput()));

        //vcf header
        conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVCFHeaderOutput() + "/vcfFileHeader.vcf");
        VariantAnnotatorEngine variantAnnotatorEngine = new VariantAnnotatorEngine(options.getAnnotationGroups(), options.getAnnotations(), null);
        VCFHeader vcfHeader = VariantCallingEngine.getVCFHeader(options, variantAnnotatorEngine, SamHdfsFileHeader.getHeader(conf));
        VCFHdfsWriter vcfHdfsWriter = new VCFHdfsWriter(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP), false, false, conf);
        vcfHdfsWriter.writeHeader(vcfHeader);
        vcfHdfsWriter.close();

        job.setJobName("GaeaGenotyper");
        job.setAnySamInputFormat(options.getInputFormat());
        conf.set(KeyIgnoringVCFOutputFormat.OUTPUT_VCF_FORMAT_PROPERTY, options.getOuptputFormat().toString());
        job.setOutputFormatClass(GaeaVCFOutputFormat.class);
        job.setOutputKeyValue(WindowsBasedWritable.class, AlignmentBasicWritable.class, NullWritable.class, VariantContextWritable.class);

        job.setJarByClass(Genotyper.class);
        if(options.getBedRegionFile() != null)
            conf.set(REFERENCE_REGION, options.getBedRegionFile());
        job.setFilterClass(GenotyperFilter.class);
        job.setWindowsBasicMapperClass(WindowsBasedAlignmentMapper.class, options.getWindowSize(),0);
        job.setReducerClass(GenotyperReducer.class);
        job.setNumReduceTasks(options.getReducerNumber());

        FileInputFormat.setInputPaths(job, new Path(options.getInput()));
        FileOutputFormat.setOutputPath(job, new Path(options.getOutput()));

        if (job.waitForCompletion(true)) {
            return 0;
        }

        return 1;
    }

    @Override
    public int run(String[] args) throws Exception {
        Genotyper genotyper = new Genotyper();
        int res = genotyper.runGenoytper(args);

        return res;
    }

}
