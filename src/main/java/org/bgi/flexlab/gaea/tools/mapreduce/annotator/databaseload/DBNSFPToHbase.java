package org.bgi.flexlab.gaea.tools.mapreduce.annotator.databaseload;

import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.HColumnDescriptor;
import org.apache.hadoop.hbase.HTableDescriptor;
import org.apache.hadoop.hbase.TableName;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.io.ImmutableBytesWritable;
import org.apache.hadoop.hbase.mapreduce.HFileOutputFormat2;
import org.apache.hadoop.hbase.mapreduce.LoadIncrementalHFiles;
import org.apache.hadoop.hbase.mapreduce.PutSortReducer;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;

import java.io.IOException;

public class DBNSFPToHbase extends ToolsRunner {
    public static final String DEFAULT_COLUMN_FAMILY = "data";
    private DBNSFPToHbaseOptions options = null;
    private Connection conn = null;
    private TableName tableName = null;

    private void createTable(TableName tableName) throws IOException {
        Admin admin = conn.getAdmin();
        if (admin.tableExists(tableName)) {
            admin.close();
            return;
        }
        HTableDescriptor table = new HTableDescriptor(tableName);
        table.addFamily(new HColumnDescriptor(DEFAULT_COLUMN_FAMILY));
        admin.createTable(table);
    }

    private void LoadHFile2HBase(Configuration conf, TableName tableName, String hfile) throws Exception {
        conf.set("hbase.metrics.showTableName", "false");
        LoadIncrementalHFiles loader = new LoadIncrementalHFiles(conf);
        Admin admin = conn.getAdmin();
        Table table = conn.getTable(tableName);
        RegionLocator rl = conn.getRegionLocator(tableName);
        loader.doBulkLoad(new Path(hfile), admin, table, rl);
    }

    @Override
    public int run(String[] args) throws Exception {
        Configuration conf = HBaseConfiguration.create();

        String[] remainArgs = remainArgs(args, conf);
        options = new DBNSFPToHbaseOptions();
        options.parse(remainArgs);
        options.setHadoopConf(remainArgs, conf);
        tableName = TableName.valueOf(options.getTableName());
        conf.set("DEFAULT_COLUMN_FAMILY", "data");
        conf.addResource(new Path(options.getConfig() + "hbase-site.xml"));
        conf.addResource(new Path(options.getConfig() + "core-site.xml"));
        conn = ConnectionFactory.createConnection(conf);

        setHeader(new Path(options.getInput()), conf);
        long reduceThreshMem = (long) (1 << 28);
        conf.setLong("putsortreducer.row.threshold", reduceThreshMem);

        Job job = Job.getInstance(conf, "dbNSFPtoHbase");
        createTable(tableName);

        job.setJarByClass(org.bgi.flexlab.gaea.tools.mapreduce.annotator.databaseload.DBNSFPToHbase.class);
        job.setMapperClass(DBNSFPToHbaseMapper.class);
        job.setReducerClass(PutSortReducer.class);

        job.setMapOutputKeyClass(ImmutableBytesWritable.class);
        job.setMapOutputValueClass(Put.class);
        job.setOutputKeyClass(ImmutableBytesWritable.class);
        job.setOutputValueClass(Put.class);

        FileInputFormat.setInputPaths(job, new Path(options.getInput()));
        FileOutputFormat.setOutputPath(job, new Path(options.getHFileOutput()));

//        HFileOutputFormat2.configureIncrementalLoad(job, new HTable(conf,options.getTableName()));
        HFileOutputFormat2.configureIncrementalLoad(job, conn.getTable(tableName), conn.getRegionLocator(tableName));

        if (job.waitForCompletion(true)) {
            LoadHFile2HBase(conf, tableName, options.getHFileOutput());
            conn.close();
            return 0;
        } else {
            conn.close();
            return 1;
        }
    }

    private void setHeader(Path inputPath, Configuration conf) {
        try {
            FileSystem fs = inputPath.getFileSystem(conf);
            fs = inputPath.getFileSystem(conf);
            if (!fs.exists(inputPath)) {
                System.out.println("Input File Path is not exist! Please check input var.");
                System.exit(-1);
            }
            if (fs.isFile(inputPath)) {
                FSDataInputStream in = fs.open(inputPath);
                AsciiLineReaderIterator it = new AsciiLineReaderIterator(new AsciiLineReader(in));
                String headerStr = it.next();
                if(!headerStr.startsWith("#")){
                    System.err.println("No header info!");
                }else {
                    String[] header = headerStr.substring(1).split("\t");
                    conf.setStrings("header", header);
                }
            } else {
                FileStatus stats[] = fs.listStatus(inputPath);
                Path filePath = stats[0].getPath();
                setHeader(filePath, conf);
            }
            fs.close();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
