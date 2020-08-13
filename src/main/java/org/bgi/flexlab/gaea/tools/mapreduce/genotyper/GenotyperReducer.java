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

import htsjdk.samtools.SAMFileHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.data.mapreduce.input.bed.RegionHdfsParser;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.AlignmentBasicWritable;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.pileup.ReadsPool;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.variant.VariantCallContext;
import org.bgi.flexlab.gaea.tools.genotyer.VariantCallingEngine;
import org.bgi.flexlab.gaea.util.Window;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.IOException;
import java.util.List;

/**
 * Created by zhangyong on 2017/3/1.
 */
public class GenotyperReducer extends Reducer<WindowsBasedWritable, AlignmentBasicWritable, NullWritable, VariantContextWritable>{
    /**
     * options
     */
    private GenotyperOptions options = new GenotyperOptions();

    /**
     * sam file header
     */
    private SAMFileHeader header;

    /**
     * shared reference
     */
    private ReferenceShare genomeShare;

    /**
     * variant calling engine
     */
    private VariantCallingEngine engine;

    /**
     * output writable
     */
    private VariantContextWritable variantContextWritable;

    /**
     * region
     */
    private RegionHdfsParser region = null;

    @Override
    protected void setup(Context context) throws IOException {
        Configuration conf = context.getConfiguration();
        options.getOptionsFromHadoopConf(conf);
        header = SamHdfsFileHeader.getHeader(conf);
        genomeShare = new ReferenceShare();
        genomeShare.loadChromosomeList(options.getReference());
        engine = new VariantCallingEngine(options, header);
        variantContextWritable = new VariantContextWritable();
        AlignmentsBasic.initIdSampleHash(header.getReadGroups());

        if (options.getBedRegionFile() != null) {
            region = new RegionHdfsParser();
            region.parseBedFileFromHDFS(options.getBedRegionFile(), false);
        }
    }

    @Override
    public void reduce(WindowsBasedWritable key, Iterable<AlignmentBasicWritable> values, Context context) throws IOException, InterruptedException {
        Window win = new Window(header, key.getChromosomeIndex(), key.getWindowsNumber(), options.getWindowSize());
        ReadsPool readsPool = new ReadsPool(values.iterator());
        engine.init(readsPool, win, genomeShare.getChromosomeInfo(header.getSequence(key.getChromosomeIndex()).getSequenceName()));

        //System.err.println("win:" + win.getContigName() + ":" + win.getStart() + "-" + win.getStop());
        List<VariantCallContext> variantContexts = engine.reduce();
        while(variantContexts != null) {
            if(variantContexts.size() == 0) {
                variantContexts = engine.reduce();
                continue;
            }
            //System.err.println("output variant before filter.");
            for (VariantCallContext vc : variantContexts) {
                if(vc.shouldEmit && vc.getStart() >= win.getStart() && vc.getStart() <= win.getStop()) {
                    if(region != null && !region.isPositionInRegion(vc.getContig(), vc.getStart() - 1)) {
                        continue;
                    }
                    //System.err.println("output variant after filter.");
                    variantContextWritable.set(vc);
                    context.write(NullWritable.get(), variantContextWritable);
                }
            }
            variantContexts = engine.reduce();
        }
        //System.err.println("finished all.");
    }

    @Override
    protected void cleanup(Context context) {

    }
}
