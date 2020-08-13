package org.bgi.flexlab.gaea.tools.mapreduce.annotator.databaseload;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hbase.client.Put;
import org.apache.hadoop.hbase.io.ImmutableBytesWritable;
import org.apache.hadoop.hbase.util.Bytes;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

import java.io.IOException;

public class DBNSFPToHbaseMapper extends Mapper<LongWritable, Text, ImmutableBytesWritable, Put> {
    private String[] headerFields;

    @Override
    protected void setup(Context context) throws IOException, InterruptedException {
        Configuration conf = context.getConfiguration();
        headerFields = conf.getStrings("header");
    }

    public void map(LongWritable key, Text value, Context context)
            throws IOException, InterruptedException {
        String vText = value.toString();
        Configuration conf = context.getConfiguration();

        String[] fields = vText.split("\t");
        if(fields[0].startsWith("#"))
            return;
//        context.getCounter("Counters", "record").increment(1);
        byte[] rowKeys = Bytes.toBytes(fields[0] + '-' + fields[1] + '-' + fields[3]);
        ImmutableBytesWritable rKey = new ImmutableBytesWritable(rowKeys);
        Put put = new Put(rowKeys);
        for (int i = 0; i < headerFields.length; i++) {
            put.addColumn(Bytes.toBytes(conf.get("DEFAULT_COLUMN_FAMILY")), Bytes.toBytes(headerFields[i]), Bytes.toBytes(fields[i]));
        }
        context.write(rKey, put);
    }
}
