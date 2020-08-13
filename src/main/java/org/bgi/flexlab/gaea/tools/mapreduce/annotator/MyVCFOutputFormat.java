package org.bgi.flexlab.gaea.tools.mapreduce.annotator;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.seqdoop.hadoop_bam.KeyIgnoringVCFOutputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.IOException;

public class MyVCFOutputFormat
        extends FileOutputFormat<Text, VariantContextWritable> {
    static final String INPUT_PATH_PROP = "vcf.input_path";

    private KeyIgnoringVCFOutputFormat<Text> baseOF;

    private void initBaseOF(Configuration conf) {
        if (baseOF == null)
            baseOF = new KeyIgnoringVCFOutputFormat<>(conf);
    }

    @Override
    public RecordWriter<Text, VariantContextWritable> getRecordWriter(
            TaskAttemptContext context)
            throws IOException {
        final Configuration conf = context.getConfiguration();
        initBaseOF(conf);

        if (baseOF.getHeader() == null) {
//            final Path p = new Path(conf.get(INPUT_PATH_PROP));
//            baseOF.readHeaderFrom(p, p.getFileSystem(conf));
        }

        return baseOF.getRecordWriter(context, getDefaultWorkFile(context, ""));
    }
}