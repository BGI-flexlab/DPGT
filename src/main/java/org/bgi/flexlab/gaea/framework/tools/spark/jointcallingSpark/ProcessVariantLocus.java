package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import scala.Tuple2;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;


public class ProcessVariantLocus implements PairFlatMapFunction<String, GenomeLongRegion, Integer> {
    private final ArrayList<GenomeLocation> regions=new ArrayList<>();
    private final Set<VCFHeaderLine> gvcfHeaderMetaInfo;
    private final String vVcfPath;
    private static final Logger logger = LoggerFactory.getLogger(ProcessVariantLocus.class);
    private VCFHeaderVersion version = null;
    private VCFEncoder vcfEncoder=null;
    private final String outputDir;
    private final DriverBC dBC;
    public static HashMap<String,BufferedReader> sampleReaders=new HashMap<>();
    public ProcessVariantLocus(ArrayList<GenomeLocation> regions, Broadcast<DriverBC> dBC) throws IOException {
        this.regions.addAll(regions);
        this.vVcfPath=dBC.value().vVcfPath;
        outputDir=dBC.value().outputDir;
        this.dBC=dBC.value();
        SeekableStream in=new SeekableFileStream(new File(vVcfPath));
        VCFHeader virtualHeader = VCFHeaderReader.readHeaderFrom(in);
        gvcfHeaderMetaInfo=virtualHeader.getMetaDataInInputOrder();
        for (final VCFHeaderLine line : gvcfHeaderMetaInfo) {
            if (VCFHeaderVersion.isFormatString(line.getKey())) {
                version = VCFHeaderVersion.toHeaderVersion(line.getValue());
                break;
            }
        }
    }

    @Override
    public Iterator<Tuple2<GenomeLongRegion,Integer>> call(String s) throws Exception {
        LinkedList<Tuple2<GenomeLongRegion,Integer>> variantsList = new LinkedList<>();
        Configuration conf=new Configuration();
        SeekableStream in2 = new SeekableFileStream(new File(outputDir+"/vcfheader"));
        VCFHeader mergedHeader = VCFHeaderReader.readHeaderFrom(in2);
        if(vcfEncoder==null) {
            vcfEncoder = new VCFEncoder(mergedHeader, true, true);
        }
        if(s.startsWith("file://"))
            s=s.substring(7);
        VCFCodec tmp_codec=new VCFCodec();
        Set<String> samples=new HashSet<>();
        samples.add(s);
        tmp_codec.setVCFHeader(new VCFHeader(gvcfHeaderMetaInfo,samples),version);

        //使用query方式试试
        VCFCodec query_codec=new VCFCodec();
        TabixFeatureReader sampleReader;

        sampleReader=new TabixFeatureReader(s,query_codec);
        logger.info("current process sample:\t"+s);
        logger.info("region size:\t"+regions.size());
        for(GenomeLocation curRegion:regions){
            logger.info(curRegion.toString());
        }
        for(GenomeLocation curGloc:regions) {
            logger.info("extract bps:\t"+curGloc.getStart()+"\t"+curGloc.getEnd());
            Iterator<VariantContext> it = sampleReader.query(curGloc.getContig(), curGloc.getStart(), curGloc.getEnd());

            //不使用query的代码结束位置
            while (it.hasNext()) {
                final VariantContext vc = it.next();
                if (vc.getNAlleles() > 2) {
                    long curPosStart=dBC.accumulateLength.get(dBC.chrIndex.get(vc.getContig()))+vc.getStart();
                    long curPosEnd=dBC.accumulateLength.get(dBC.chrIndex.get(vc.getContig()))+vc.getEnd();
                    variantsList.add(new Tuple2<>(new GenomeLongRegion(curPosStart,curPosEnd), 1));
                }
            }

        }
        return variantsList.iterator();
    }

}
