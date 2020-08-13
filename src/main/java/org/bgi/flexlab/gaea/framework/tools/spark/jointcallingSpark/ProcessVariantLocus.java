package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import scala.Tuple2;

import java.io.*;
import java.util.*;


public class ProcessVariantLocus implements PairFlatMapFunction<String, GenomeLongRegion, Integer> {
    private final GenomeLongRegion processRegion;
    private final ArrayList<GenomeLocation> regions=new ArrayList<>();
    private Set<VCFHeaderLine> gvcfHeaderMetaInfo;
    private final String vVcfPath;
    private static final Logger logger = LoggerFactory.getLogger(ProcessVariantLocus.class);
    private VCFHeaderVersion version = null;
    private VCFEncoder vcfEncoder=null;
    private VCFHeader mergedHeader=null;
    private final String outputDir;
    private final DriverBC dBC;
    public static HashMap<String,BufferedReader> sampleReaders=new HashMap<>();
    public ProcessVariantLocus(GenomeLongRegion region,ArrayList<GenomeLocation> regions,Broadcast<DriverBC> dBC) throws IOException {
        this.processRegion=region;
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
//        VCFFileReader fileReader=dBC.sampleReader.get(s);
//        Iterator<VariantContext> vcsAtRegion=fileReader.query(processRegion.getContig(),processRegion.getStart(),processRegion.getEnd());
        Configuration conf=new Configuration();
        SeekableStream in2 = new SeekableFileStream(new File(outputDir+"/vcfheader"));
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
        //先处理lastVC
        //若使用query，则不需要处理lastVC
        //使用lastVC的代码起始位置
        /*
        VariantContext lastVc=null;
        File lastVCFile=new File((new File(vVcfPath)).getParent()+"/"+s+".mr1.lastVC");
        if(lastVCFile.exists()){
            BufferedReader lastVCReader=new BufferedReader(new FileReader(lastVCFile));
            String line1Sample,line2VC;
            if((line1Sample=lastVCReader.readLine())!=null){
                if(!line1Sample.equals(s)){
                    logger.error("code error");
                    System.exit(1);
                }
            }
            if((line2VC=lastVCReader.readLine())!=null){
                lastVc=tmp_codec.decode(line2VC);
            }
            lastVCReader.close();
        }
        if(lastVc!=null){
            if(processRegion.getContigIndex()==dBC.chrIndex.get(lastVc.getContig())){
                if(processRegion.getStart()<=lastVc.getEnd() && processRegion.getEnd()>=lastVc.getStart()){
                    GenomeLongRegion gloc=new GenomeLongRegion(lastVc.getContig(),dBC.chrIndex.get(lastVc.getContig()),lastVc.getStart(),lastVc.getEnd());
                    variantsList.add(new Tuple2<>(gloc,1));
                }
            }
        }
        //使用lastVC的代码结束位置
        */
        //使用query方式试试
        File fPath=new File(s);
        VCFCodec query_codec=new VCFCodec();
        Path idxFile=new Path(s+".tbi");
        FileSystem idxFs=idxFile.getFileSystem(conf);
        TabixFeatureReader sampleReader=null;
//        if(!idxFs.exists(idxFile)){
//            TabixIndex idx= IndexFactory.createTabixIndex(new File(s),query_codec,VCF,mergedHeader.getSequenceDictionary());
//            LittleEndianOutputStream stream = new LittleEndianOutputStream(new BufferedOutputStream(new BlockCompressedOutputStream(new File(s+".tbi"))));
//            idx.write(stream);
//            stream.close();
//        }
        sampleReader=new TabixFeatureReader(s,query_codec);
        logger.info("current process sample:\t"+s);
        logger.info("region size:\t"+regions.size());
        for(GenomeLocation curRegion:regions){
            logger.info(curRegion.toString());
        }
        for(GenomeLocation curGloc:regions) {
            logger.info("extract bps:\t"+curGloc.getStart()+"\t"+curGloc.getEnd());
            Iterator<VariantContext> it = sampleReader.query(curGloc.getContig(), curGloc.getStart(), curGloc.getEnd());
            //不使用query的代码起始位置
//        Integer sampleIdx=dBC.sampleIndex.get(dBC.pathSample.get(fPath.getName()));
//        while((vcfLine=JointCallingSpark.sampleReaders.get(sampleIdx).readLine())!=null){
            //不使用query的代码结束位置
            while (it.hasNext()) {
                final VariantContext vc = it.next();
//            if(vcfLine.startsWith("#")) {
//                continue;
//            }
//            final VariantContext vc = tmp_codec.decode(vcfLine);
                if (vc.getNAlleles() > 2) {
                    long curPosStart=dBC.accumulateLength.get(dBC.chrIndex.get(vc.getContig()))+vc.getStart();
                    long curPosEnd=dBC.accumulateLength.get(dBC.chrIndex.get(vc.getContig()))+vc.getEnd();
                    variantsList.add(new Tuple2<>(new GenomeLongRegion(curPosStart,curPosEnd), 1));
//                    if (curGloc.getContigIndex() < dBC.chrIndex.get(vc.getContig())) {
//                        //writeToDisk(s,vc); 若使用query方式，这里不需要保存
//                        break;
//                    } else if (curGloc.getContigIndex() > dBC.chrIndex.get(vc.getContig())) {
//                        continue;
//                    } else {
//                        if (curGloc.getStart() > vc.getEnd()) {
//                            continue;
//                        } else if (curGloc.getEnd() < vc.getStart()) {
//                            //因为不能回滚，无法保存当前值，而下个处理区间还需要用，所以要保存起来
//                            //累加器不可靠，这里只能借助外存储了，将其写到临时文件里面
//                            //写的格式为（两行）:文件名（其实是路径）\nVariantContext.toString()
//                            //writeToDisk(s,vc); //若使用query方式，这里不需要保存
//                            break;
//                        } else {
//                            GenomeLongRegion gloc = new GenomeLongRegion(vc.getContig(), dBC.chrIndex.get(vc.getContig()), vc.getStart(), vc.getEnd());
//                            variantsList.add(new Tuple2<>(gloc, 1));
//                        }
//                    }
                }
            }
//        while(vcsAtRegion.hasNext()){
//            VariantContext vc=vcsAtRegion.next();
//            if(vc.getNAlleles()>2) {
//                GenomeLongRegion gloc=new GenomeLongRegion(vc.getContig(),dBC.chrIndex.get(vc.getContig()),vc.getStart(),vc.getEnd());
//                variantsList.add(new Tuple2<>(gloc,1));
//            }
//        }
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
