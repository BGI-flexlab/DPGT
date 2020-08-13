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

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;
import org.bgi.flexlab.gaea.tools.vcfstats.report.PerSampleVCFReport;
import org.bgi.flexlab.gaea.tools.vcfstats.report.VCFReport;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.IOException;

public class VCFStatsMapper extends Mapper<LongWritable, VariantContextWritable, NullWritable, Text> {

//    private DbsnpShare dbsnpShare = null;

    /**
     * vcf loader
     */
//    private VCFLocalLoader DBloader = null;

    private VCFReport vcfReport;
    private Text resultValue = new Text();
    protected void setup(Context context)
            throws IOException, InterruptedException {
        VCFStatsOptions options = new VCFStatsOptions();
        options.getOptionsFromHadoopConf(context.getConfiguration());
        vcfReport = new VCFReport(options);

//        if(options.getDbsnpFile() != null) {
//            dbsnpShare = new DbsnpShare(options.getDbsnpFile(), options.getReferenceSequencePath());
//            dbsnpShare.loadChromosomeList(options.getDbsnpFile() + VcfIndex.INDEX_SUFFIX);
//            DBloader = new VCFLocalLoader(options.getDbsnpFile());
//        }
    }

    private boolean filteVariant(VariantContext variantContext){
        return !variantContext.isVariant();
    }

    @Override
    protected void map(LongWritable key, VariantContextWritable value, Context context) throws IOException, InterruptedException {
        VariantContext vc = value.get();
        if(filteVariant(vc))
            return;
        vcfReport.parseVariation(vc);
    }

    @Override
    protected void cleanup(Context context) throws IOException, InterruptedException {
        for(String sample: vcfReport.getPerSampleVCFReports().keySet()){
            PerSampleVCFReport sampleVCFReport = vcfReport.getPerSampleVCFReports().get(sample);
            resultValue.set(sampleVCFReport.toReducerString());
            context.write(NullWritable.get(), resultValue);
        }
    }
}
