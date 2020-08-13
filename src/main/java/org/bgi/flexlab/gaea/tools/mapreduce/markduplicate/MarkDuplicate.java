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
 * Copyright (c) 2009-2012 The Broad Institute
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
package org.bgi.flexlab.gaea.tools.mapreduce.markduplicate;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.bam.GaeaBamInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.bam.GaeaBamOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.DuplicationKeyWritable;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.seqdoop.hadoop_bam.SAMFormat;

import java.io.IOException;

public class MarkDuplicate extends ToolsRunner{

    public MarkDuplicate(){
        this.toolsDescription = "Gaea Mark PCR duplication";
    }

    public int runMarkDuplicate(String[] args) throws IOException, ClassNotFoundException, InterruptedException {
        BioJob job = BioJob.getInstance();
        Configuration conf = job.getConfiguration();
        String[] remainArgs = remainArgs(args, conf);
        MarkDuplicateOptions options = new MarkDuplicateOptions();
        options.parse(remainArgs);
        options.setHadoopConf(remainArgs, conf);

        job.setHeader(new Path(options.getInput()), new Path(options.getOutput()));

        job.setJobName("GaeaMarkDuplicate");
        job.setJarByClass(MarkDuplicate.class);
        job.setMapperClass(MarkDuplicateMappper.class);
        job.setReducerClass(MarkDuplicateReducer.class);
        job.setNumReduceTasks(options.getReducerNum());

        job.setMapOutputKeyClass(DuplicationKeyWritable.class);
        job.setMapOutputValueClass(SamRecordWritable.class);

        job.setOutputKeyClass(NullWritable.class);
        job.setOutputValueClass(SamRecordWritable.class);

        job.setAnySamInputFormat(options.getInputFormat());
        if(options.getOutputFormat() == 0){
            job.setOutputFormatClass(GaeaBamOutputFormat.class);
        }

        FileInputFormat.setInputPaths(job, options.getInputFileList().toArray(new Path[options.getInputFileList().size()]));
        Path oPath = new Path(options.getOutput()+"/Mark");
        FileOutputFormat.setOutputPath(job, oPath);
        boolean success = job.waitForCompletion(true);
        return success ? 0 : 1;
    }

    @Override
    public int run(String[] args) throws Exception {
        MarkDuplicate md = new MarkDuplicate();
        return md.runMarkDuplicate(args);
    }
}
