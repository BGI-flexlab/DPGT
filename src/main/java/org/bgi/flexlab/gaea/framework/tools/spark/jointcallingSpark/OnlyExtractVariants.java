package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.util.SerializableConfiguration;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import scala.Tuple2;

import java.io.*;
import java.net.URI;
import java.util.*;
import java.util.zip.GZIPInputStream;


public class OnlyExtractVariants implements PairFlatMapFunction<String, GenomeLocation, Integer> {
    public GenomeLocation processRegion;
    public Set<VCFHeaderLine> gvcfHeaderMetaInfo;
    public String vVcfPath;
    public static Logger logger = LoggerFactory.getLogger(ProcessVariantLocus.class);
    public VCFHeaderVersion version = null;
    public VCFEncoder vcfEncoder=null;
    public VCFHeader mergedHeader=null;
    public String outputDir;
    public DriverBC dBC;
    public static HashMap<String,BufferedReader> sampleReaders=new HashMap<>();
    public OnlyExtractVariants(GenomeLocation region, Broadcast<DriverBC> dBC) throws IOException {
        this.processRegion=region;
        this.vVcfPath=dBC.value().vVcfPath;
        outputDir=dBC.value().outputDir;
        this.dBC=dBC.value();
        SeekableStream in=new SeekableFileStream(new File(vVcfPath));
        VCFHeader virtualHeader = VCFHeaderReader.readHeaderFrom(in);
        gvcfHeaderMetaInfo=virtualHeader.getMetaDataInInputOrder();
        if(version==null) {
            for (final VCFHeaderLine line : gvcfHeaderMetaInfo) {
                if (VCFHeaderVersion.isFormatString(line.getKey())) {
                    version = VCFHeaderVersion.toHeaderVersion(line.getValue());
                    break;
                }
            }
        }
    }

    @Override
    public Iterator<Tuple2<GenomeLocation,Integer>> call(String s) throws Exception {
        if(!sampleReaders.containsKey(s)) {
            Configuration hadoopConf = new Configuration();
            BufferedReader gvcfReader = null;
            if (s.startsWith("file://")) {
                s = s.substring(7);
                if (s.endsWith(".gz")) {
                    gvcfReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(s))));
                } else {
                    gvcfReader = new BufferedReader(new FileReader(s));
                }
            } else {
                FileSystem fs = FileSystem.get(URI.create(s), hadoopConf);
                FSDataInputStream fsr = fs.open(new Path(s));
                if (s.endsWith(".gz")) {
                    gvcfReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(fsr)));
                } else {
                    gvcfReader = new BufferedReader(new InputStreamReader(fsr));
                }
            }
        }
        List<Tuple2<GenomeLocation,Integer>> variantsList = new ArrayList<>();
        Configuration conf=new Configuration();
        Path path = new Path(outputDir+"/vcfheader");
        SeekableStream in2 = WrapSeekable.openPath(path.getFileSystem(conf), path);
        mergedHeader = VCFHeaderReader.readHeaderFrom(in2);
        if(vcfEncoder==null) {
            vcfEncoder = new VCFEncoder(mergedHeader, true, true);
        }
        if(s.startsWith("file://"))
            s=s.substring(7);
        String vcfLine;
        VCFCodec tmp_codec=new VCFCodec();
        Set<String> samples=new HashSet<>();
        samples.add(s);
        tmp_codec.setVCFHeader(new VCFHeader(gvcfHeaderMetaInfo,samples),version);
        File fPath=new File(s);
        VCFCodec query_codec=new VCFCodec();
        TabixFeatureReader sampleReader=new TabixFeatureReader(s,query_codec);
        logger.info("current process sample:\t"+s);
        Iterator<VariantContext> it=sampleReader.query(processRegion.getContig(),processRegion.getStart(),processRegion.getEnd());
        while(it.hasNext()){
            final VariantContext vc=it.next();
            GenomeLocation gloc=new GenomeLocation(vc.getContig(),dBC.chrIndex.get(vc.getContig()),vc.getStart(),vc.getEnd());
            variantsList.add(new Tuple2<>(gloc,1));
        }

        return variantsList.iterator();
    }

    private void writeToDisk(String s, VariantContext vc) throws IOException {
        //写到哪里
        FileWriter w=new FileWriter((new File(vVcfPath)).getParent()+"/"+s+".mr1.lastVC");
        w.write(s);
        w.write("\n");
        w.write(vcfEncoder.encode(vc));
        w.write("\n");
        w.close();
        //最后在driver中要清理
    }
}
