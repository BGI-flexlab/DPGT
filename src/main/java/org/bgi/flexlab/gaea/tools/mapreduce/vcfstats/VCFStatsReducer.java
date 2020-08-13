package org.bgi.flexlab.gaea.tools.mapreduce.vcfstats;

import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.output.MultipleOutputs;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.IOException;

public class VCFStatsReducer extends Reducer<LongWritable,VariantContextWritable,
            NullWritable,VariantContextWritable> {
    private MultipleOutputs<NullWritable, Text> mos;

    @Override
    protected void setup(Context context) throws IOException, InterruptedException {
        super.setup(context);
    }

    @Override
    protected void reduce(LongWritable ignored, Iterable<VariantContextWritable> records, Context ctx)
            throws IOException, InterruptedException {
        for (VariantContextWritable rec : records)
            ctx.write(NullWritable.get(), rec);
    }

    @Override
    protected void cleanup(Context context) throws IOException, InterruptedException {
        super.cleanup(context);
    }
}
