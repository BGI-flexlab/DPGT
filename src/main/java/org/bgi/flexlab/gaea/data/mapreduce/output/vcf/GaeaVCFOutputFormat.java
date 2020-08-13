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
package org.bgi.flexlab.gaea.data.mapreduce.output.vcf;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.seqdoop.hadoop_bam.KeyIgnoringVCFOutputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.IOException;

/**
 * Created by zhangyong on 2017/3/3.
 * came form VCFsort
 */
public class GaeaVCFOutputFormat<K> extends FileOutputFormat<K, VariantContextWritable> {

    public static final String OUT_PATH_PROP = "gaea.vcf.outpath";
    public static final String HEADER_MODIFY = "gaea.vcf.header.modify";

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
        	if(conf.get(OUT_PATH_PROP) != null){
        		final Path p = new Path(conf.get(OUT_PATH_PROP));
        		baseOF.readHeaderFrom(p, p.getFileSystem(conf));
        	}
        }
        
        if(conf.getBoolean(GaeaVCFOutputFormat.HEADER_MODIFY, false)){
        	final boolean wh = context.getConfiguration().getBoolean(
        			KeyIgnoringVCFOutputFormat.WRITE_HEADER_PROPERTY, true);
        	return new GaeaKeyIgnoringVCFRecordWriter<K>(getDefaultWorkFile(context, ""),baseOF.getHeader(),wh,context);
        }

        return baseOF.getRecordWriter(context, getDefaultWorkFile(context, ""));
    }

    // Allow the output directory to exist.
    @Override public void checkOutputSpecs(JobContext job) {}
}
