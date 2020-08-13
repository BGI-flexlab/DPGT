package org.bgi.flexlab.gaea.tools.mapreduce.fastqmerge;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.hbase.util.Bytes;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.CompressionOutputStream;
import org.apache.hadoop.util.ReflectionUtils;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;

import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 * Load Clean Fastq Data from HDFS!
 *
 */
public class FastqMerge  extends ToolsRunner {

    private Configuration conf;
    private FastqMergeOptions options;

    public FastqMerge() {
        this.toolsDescription = "Load Clean Fastq Data from HDFS";
    }

    private int runFastqMerge(String[] arg0) throws Exception {

        conf = new Configuration();
        String[] remainArgs = remainArgs(arg0, conf);

        options = new FastqMergeOptions();
        options.parse(remainArgs);
        options.setHadoopConf(remainArgs, conf);

        Path rlist = new Path(options.getInput());
        Path fq1Path = new Path(options.getFq1out());
        FileSystem outFS = fq1Path.getFileSystem(conf);
        FSDataOutputStream fq1os = outFS.create(fq1Path);

        FSDataOutputStream fq2os = null;

        if(options.getFq2out() != null){
            Path fq2Path = new Path(options.getFq2out());
            fq2os = outFS.create(fq2Path);
        }

        FileSystem fs = FileSystem.get(conf);

        Class<?> codecClass = Class.forName("org.apache.hadoop.io.compress.GzipCodec");

        CompressionCodec codec = (CompressionCodec) ReflectionUtils.newInstance(codecClass, conf);

        CompressionOutputStream fq1GZos = null;
        CompressionOutputStream fq2GZos = null;

        boolean compress = false;

        if(options.getFq1out().endsWith(".gz")){
            compress = true;
            fq1GZos = codec.createOutputStream(fq1os);
            fq2GZos = codec.createOutputStream(fq2os);
        }

        Path[] listedPaths = FileUtil.stat2Paths(fs.listStatus(rlist));
        for (Path p : listedPaths) {

            String[] basename = p.toString().split("/");

            if (!basename[basename.length-1].startsWith("part-")) {
                continue;
            }

            System.out.println(p);

            FSDataInputStream dis = fs.open(p);
            BufferedReader br = new BufferedReader(new InputStreamReader(dis));


            String line;
            while((line = br.readLine()) != null){
                String bqHead, fastqStr;

                line = line.trim();
                fastqStr = line + "\n";
                fastqStr += br.readLine().trim() + "\n";
                bqHead = br.readLine();
                bqHead = "+";
                fastqStr +=  bqHead + "\n";
                fastqStr += br.readLine().trim() + "\n";

                if (line.endsWith("/1")) {
                    if(compress)
                        fq1GZos.write(Bytes.toBytes(fastqStr));
                    else
                        fq1os.write(Bytes.toBytes(fastqStr));
                }
                else if (line.endsWith("/2")) {
                    assert options.getFq2out() != null : "The Fastq is PE data, Please set the -D,--cleanFq2 parameter!";
                    if(compress)
                        fq2GZos.write(Bytes.toBytes(fastqStr));
                    else
                        fq2os.write(Bytes.toBytes(fastqStr));
                }else {
                    assert options.getFq2out() == null;
                    if(compress)
                        fq1GZos.write(Bytes.toBytes(fastqStr));
                    else
                        fq1os.write(Bytes.toBytes(fastqStr));
                }
            }
        }

        if(compress) {
            fq1GZos.close();
            fq2GZos.close();
        }
        else {
            fq1os.close();
            if(options.getFq2out() != null) {
                fq2os.close();
            }
        }

        return 1;
    }

    @Override
    public int run(String[] args) throws Exception {
        FastqMerge fastqMerge = new FastqMerge();
        return fastqMerge.runFastqMerge(args);
    }
}
