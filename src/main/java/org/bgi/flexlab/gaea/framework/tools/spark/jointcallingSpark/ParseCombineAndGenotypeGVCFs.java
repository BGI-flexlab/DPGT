
package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.util.LongAccumulator;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.structure.dbsnp.DbsnpShare;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalWriter;
import org.bgi.flexlab.gaea.data.variant.filter.VariantRegionFilter;
import org.bgi.flexlab.gaea.tools.jointcalling.JointCallingEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;

public class ParseCombineAndGenotypeGVCFs implements Function2<Integer,Iterator<String>,Iterator<String>> {
    private HashMap<Integer, String> contigs = null;
    private JointCallingEngine engine = null;
//    private HashMap<String,Integer> chrIndexs = new HashMap<>();
    private GenomeLocationParser parser = null;
    private ReferenceShare genomeShare = null;
    private DbsnpShare dbsnpShare = null;
    private VCFLocalLoader loader = null;
    private VariantRegionFilter filter = null;
    private VCFHeader header = null;
    private VCFHeader header2=null;
    private final MultipleVCFHeaderForJointCalling headers = new MultipleVCFHeaderForJointCalling();
//    private ArrayList<ArrayList<String>> multiMapSampleNames=new ArrayList<ArrayList<String> >();
    private VCFEncoder vcfEncoder=null;
//    public BufferedReader bp_reader=null;
private String winLine=null;
    private GenomeLongRegion region=null;
    private final SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
    private HashMap<String, String> confMap=new HashMap<>();
    private final String bpFile;
    private final String[] args;
    private final DriverBC dBC;
    private ArrayList<GenomeLocation> regions=new ArrayList<>();
    private final int cycleIter;
    private final ArrayList<Long> bpPartition=new ArrayList<>();
    private final LongAccumulator totalVariantsNum;
    public ParseCombineAndGenotypeGVCFs(GenomeLongRegion region, ArrayList<GenomeLocation> regions, String[] args, String outputBP, HashMap<String, String> confMap,
                                        Broadcast<DriverBC> dBC, int cycleIter, ArrayList<Long> bpPartition,
                                        LongAccumulator totalVariantsNum) {
        this.args=args;
        this.dBC=dBC.value();
        this.cycleIter=cycleIter;
        this.confMap=confMap;
        this.region=region;
        this.regions=regions;
        bpFile=outputBP;
        this.totalVariantsNum=totalVariantsNum;
        this.bpPartition.addAll(bpPartition);
    }
    class VcComp implements Comparator<VariantContext>,Serializable{

        @Override public int compare(VariantContext o1, VariantContext o2) {
            return o1.getStart()-o2.getStart();
        }
    }
    private final Comparator<VariantContext> comparator4 = new VcComp();
    String HEADER_DEFAULT_PATH = "vcfheader";
    private final String MERGER_HEADER_INFO = "vcfheaderinfo";
    @Override public Iterator<String> call(Integer index,Iterator<String> stringIterator) throws IOException {
        Logger log= LoggerFactory.getLogger(ParseCombineAndGenotypeGVCFs.class);
        Configuration hadoop_conf=new Configuration();
        for(String k:confMap.keySet()){
            hadoop_conf.set(k,confMap.get(k));
        }
//        ArrayList<String> mapGvcfList=new ArrayList<String>();

//        while(stringIterator.hasNext()) {
//            String gvcfPath = stringIterator.next();
//            mapGvcfList.add(gvcfPath);
//        }
        contigs = new HashMap<>();
//        SeekableStream in=new SeekableFileStream(new File(confMap.get(GaeaVCFOutputFormat.OUT_PATH_PROP)));
//		SeekableStream in=new SeekableFileStream(new File(options.getOutDir()()+"/virtual.vcf"));
        System.out.println(formatter.format(new Date())+"\tbefore readHeader");
//        Path path = new Path(confMap.get(GaeaVCFOutputFormat.OUT_PATH_PROP));
        SeekableStream in = new SeekableFileStream(new File(dBC.outputDir+"/vcfheader"));
        header = VCFHeaderReader.readHeaderFrom(in);
        System.out.println(formatter.format(new Date())+"\tafter readHeader");
        in.close();
        if(header == null)
            throw new RuntimeException("header is null !!!");

        for (VCFContigHeaderLine line : header.getContigLines()) {
            contigs.put(line.getContigIndex(), line.getID());
//            chrIndexs.put(line.getID(),line.getContigIndex());
        }
        parser = new GenomeLocationParser(header.getSequenceDictionary());
        hadoop_conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, dBC.options.getOutDir() + "/vcfheader");
        hadoop_conf.set(MERGER_HEADER_INFO, dBC.options.getOutDir()+"/"+MERGER_HEADER_INFO);
//        headers.readHeaders(hadoop_conf);
        SeekableStream in2=new SeekableFileStream(new File(dBC.outputDir+"/vcfheader"));
        headers.setMergeHeader(VCFHeaderReader.readHeaderFrom(in2));
        headers.setCurrentIndex(dBC.sampleIndex.size());
        HashMap<Integer, Set<String>> tmpNameHeaders=new HashMap<Integer, Set<String>>();
        for(Map.Entry<String,Integer> kv:dBC.sampleIndex.entrySet()){
            TreeSet<String> samples=new TreeSet<>();
            samples.add(kv.getKey());
            tmpNameHeaders.put(kv.getValue(),samples);
        }
        headers.setNameHeaders(tmpNameHeaders);
        String sampleStr = confMap.get(dBC.INPUT_ORDER);
        String[] allSample=sampleStr.split(",");
//        int mapperLine=allSample.length/options.getMapperNumber()<2?2:allSample.length/options.getMapperNumber();
//        int mapperLine=10;


        System.out.println(formatter.format(new Date())+"\tbefore engine init");
        dBC.options.setS(dBC.options.STANDARD_CONFIDENCE_FOR_CALLING);
        dBC.options.sets(dBC.options.STANDARD_CONFIDENCE_FOR_EMITTING);
        engine = new JointCallingEngine(dBC.options, parser,header,headers,allSample,dBC.multiMapSampleNames,hadoop_conf);
        System.out.println(formatter.format(new Date())+"\tafter engine init");
        genomeShare = new ReferenceShare();
        genomeShare.loadChromosomeList(dBC.options.getReference());
        dbsnpShare = new DbsnpShare(dBC.options.getDBSnp(), dBC.options.getReference());
        dbsnpShare.loadChromosomeList(dBC.options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
        loader = new VCFLocalLoader(dBC.options.getDBSnp());
        filter = new VariantRegionFilter();
        header2=engine.getVCFHeader();
        if(header2 == null)
            throw new RuntimeException("header is null !!!");
//        Path pBPpath=new Path(bpFile);
        String mergedHeaderFile=dBC.options.getOutDir()+"/merged.vcfheader";
        File hFile=new File(mergedHeaderFile);
        if(!hFile.exists()) {
            VCFLocalWriter mergedVCFHeaderWriter = new VCFLocalWriter(mergedHeaderFile, false, false);
            mergedVCFHeaderWriter.writeHeader(header2);
            mergedVCFHeaderWriter.close();
        }
        BufferedReader bp_reader=new BufferedReader(new FileReader(bpFile));
        winLine=bp_reader.readLine();
        System.out.println("reduce setup done");

        System.out.println(formatter.format(new Date())+"\treduce start");
//        ArrayList<BufferedReader> samplesReader=new ArrayList<>();
        ArrayList<BlockCompressedInputStream> samplesReader=new ArrayList<>();
        ArrayList<VCFCodec> samplesCodec=new ArrayList<>();
        ArrayList<Integer> samplesTag=new ArrayList<>();
//        ArrayList<Long> startOffset=new ArrayList<>();
        ArrayList<Long> endOffset=new ArrayList<>();
        for(ArrayList<String> samplesInMap:dBC.multiMapSampleNames){
            //打开文件句柄
            Integer mapSMtagInt=0;
            int samplesInMapIndex=0;
            for(String sampleName:samplesInMap){
                if(samplesInMapIndex==0) {
                    mapSMtagInt = dBC.sampleIndex.get(sampleName);
                }
                samplesInMapIndex++;
            }
            mapSMtagInt+=headers.getHeaderSize();
            samplesTag.add(mapSMtagInt);
//            String targetFileName=dBC.options.getOutDir()+"/combine."+cycleIter+"/combineFile."+mapSMtagInt+"."+index;
            String targetFileName=dBC.options.getOutDir()+"/combine."+cycleIter+"/combineFile."+mapSMtagInt;
            System.out.println("processed file:\t"+targetFileName);
//            Path filePath=new Path(targetFileName);
//            FileSystem fileFs=filePath.getFileSystem(hadoop_conf);
//            if(!fileFs.exists(filePath)){
//                log.error("code error1\t"+filePath.getName()+" not exists");
//                System.exit(1);
//            }

//            BufferedReader reader=new BufferedReader(new InputStreamReader(fileFs.open(filePath)));
//            BufferedReader reader=new BufferedReader(new FileReader(targetFileName));
            BlockCompressedInputStream reader=new BlockCompressedInputStream(new File(targetFileName));
            BufferedReader idxReader=new BufferedReader(new FileReader(targetFileName+".idx"));
            String idxLine;
            boolean findIdx=false;
            while((idxLine=idxReader.readLine())!=null){
                String[] eles=idxLine.split("\t");
                int idx=Integer.parseInt(eles[0]);
                if(idx==index) {
                    reader.seek(Long.parseLong(eles[1]));
                    endOffset.add(Long.parseLong(eles[2]));
                    findIdx=true;
                    break;
                }
            }
            idxReader.close();
            if(!findIdx){
                reader.seek(0);
                endOffset.add(-1L);
            }

            VCFHeader curSamplesMergedHeader=new VCFHeader(dBC.virtualHeader.getMetaDataInInputOrder(),samplesInMap);
            VCFCodec codec=new VCFCodec();
            codec.setVCFHeader(curSamplesMergedHeader,dBC.version);
            samplesCodec.add(codec);
            samplesReader.add(reader);
        }
        if(samplesReader.size()!=dBC.multiMapSampleNames.size() || samplesCodec.size()!=dBC.multiMapSampleNames.size()){
            log.error("code error2\t"+samplesReader.size()+"\t"+dBC.multiMapSampleNames.size());
            System.exit(1);
        }
        int vcNum=0;
        String outFile=dBC.options.getOutDir()+"/genotype."+cycleIter+"/"+region.getStart()+"_"+region.getEnd()+"."+index+".gz";
        Path outPath=new Path(outFile);
        FileSystem outFs=outPath.getFileSystem(hadoop_conf);
        if(outFs.exists(outPath)){
            outFs.delete(outPath,true);
        }
        BlockCompressedOutputStream out=new BlockCompressedOutputStream(outFile);
        long rangeInEachOutFile=(region.getEnd()-region.getStart())/dBC.options.getReducerNumber();
        log.info("range in each out file:\t"+rangeInEachOutFile);
        log.info("current index:\t"+index);
//        int minRange=100000;
//        if(rangeInEachOutFile<minRange){
//            rangeInEachOutFile=minRange;
//        }
        boolean[] readerDone=new boolean[dBC.multiMapSampleNames.size()];
        for(int i=0;i<dBC.multiMapSampleNames.size();i++){
            readerDone[i]=false;
        }
        vcfEncoder = new VCFEncoder(header, true, true);
        for(GenomeLocation curRegion:regions) {
            int start = curRegion.getStart();
            String chr = curRegion.getContig();
            int chrInx = dBC.chrIndex.get(chr);
//            if(chrInx>=25){
//                continue;
//            }
            int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
            int end = curRegion.getEnd();
            long regionStart=dBC.accumulateLength.get(chrInx)+start;
            long regionEnd=dBC.accumulateLength.get(chrInx)+end;
            if(regionStart>bpPartition.get(index)){
                break;
            }else{
                if(index>0){
                    if(regionEnd<bpPartition.get(index-1)){
                        continue;
                    }
                }
            }
            ArrayList<VariantContext> dbsnps = new ArrayList<>();
//            int realStart=0;
//            if(index==0){
//                realStart=start;
//            }else if(index>0){
//                realStart= (int) (bpPartition.get(index-1)+1-dBC.accumulateLength.get(chrInx));
//            }
//
//            int dbSnpRegionEnd=(int) (bpPartition.get(index)-dBC.accumulateLength.get(chrInx));
//
//            if(dbSnpRegionEnd>end){
//                dbSnpRegionEnd= end;
//            }
//            long startPosition = dbsnpShare.getStartPosition(chr,realStart / dBC.options.getWindowsSize(),dbSnpRegionEnd/dBC.options.getWindowsSize(),dBC.options.getWindowsSize());
//            if (startPosition >= 0) dbsnps = filter.loadFilter(loader, chr, startPosition, dbSnpRegionEnd);
//            engine.init(dbsnps);

//            if(dbsnps.size()>0) {
//                System.out.println(dbsnps.get(0));
//                System.out.println(dbsnps.get(dbsnps.size() - 1));
//            }
            //为了避免每次从头读取winBp，保留bp_reader，使其只读一遍
            System.out.println(formatter.format(new Date()) + "\tbp window before get bps:\t" + winLine);
            LinkedList<VariantContext> regionVcs=new LinkedList<>();
//            long smallWinStartLong=dBC.accumulateLength.get(chrInx)+start+rangeInEachOutFile*index;
            long smallWinStartLong=regionStart;
            long smallWinEndLong=smallWinStartLong+dBC.options.getWindowsSize()-1;
            if(smallWinEndLong>regionEnd){
                smallWinEndLong=regionEnd;
            }
            int smallWinStart=start;
            int smallWinEnd=smallWinStart+dBC.options.getWindowsSize()-1;
            if(smallWinEnd>contigLength){
                smallWinEnd=contigLength;
            }
            if(smallWinEndLong-smallWinStartLong>contigLength || smallWinEnd-smallWinStart>contigLength){
                log.error("code error");
                System.exit(1);
            }
            HashMap<Integer,ArrayList<VariantContext>> lastVCs=new HashMap<>();
//            ArrayList<BufferedWriter> outWriterList=new ArrayList<>();
            int logIter=0;
            while(smallWinStartLong<=regionEnd){
                //一个个小窗口开始处理
                //获取小窗口内的所有变异
                if(smallWinStartLong>bpPartition.get(index)){
                    break;
                }
                if(index>0){
                    if(smallWinEndLong<bpPartition.get(index-1)) {
                        smallWinStartLong=smallWinEndLong+1;
                        smallWinEndLong=smallWinStartLong+dBC.options.getWindowsSize()-1;
                        if(smallWinEndLong>regionEnd){
                            smallWinEndLong=regionEnd;
                        }
                        smallWinStart=smallWinEnd+1;
                        smallWinEnd=smallWinStart+dBC.options.getWindowsSize()-1;
                        if(smallWinEnd>contigLength){
                            smallWinEnd=contigLength;
                        }
                        if(smallWinEndLong-smallWinStartLong>contigLength || smallWinEnd-smallWinStart>contigLength){
                            log.error("code error");
                            System.exit(1);
                        }
                        continue;
                    }
                }else{
                    if(smallWinEndLong>=bpPartition.get(index)){
                        smallWinEndLong=bpPartition.get(index)-1;
                    }
                }
                if(smallWinStartLong>smallWinEndLong){
                    break;
                }
                long startPosition = dbsnpShare.getStartPosition(chr,smallWinStart / dBC.options.getWindowsSize(),smallWinEnd/dBC.options.getWindowsSize(),dBC.options.getWindowsSize());
                if (startPosition >= 0) dbsnps = filter.loadFilter(loader, chr, startPosition, smallWinEnd);
                engine.init(dbsnps);
                System.out.println("dbsnp size:\t"+dbsnps.size());
                Set<Integer> bps=new TreeSet<>();
                for(int i=0;i<dBC.multiMapSampleNames.size();i++){
                    if(readerDone[i]){
                        continue;
                    }
                    boolean remove=false;
                    boolean breakFor=false;
                    if(lastVCs.containsKey(i)){
                        for(VariantContext lastVC:lastVCs.get(i)){
                            long curPosStartLong=dBC.accumulateLength.get(dBC.chrIndex.get(lastVC.getContig()))+lastVC.getStart();
                            long curPosEndLong=dBC.accumulateLength.get(dBC.chrIndex.get(lastVC.getContig()))+lastVC.getEnd();
                            if(curPosStartLong<=smallWinEndLong && curPosEndLong>=smallWinStartLong){
                                remove=true;
                                if(lastVC.getNAlleles()>2) {
                                    for (int pos = lastVC.getStart(); pos <= lastVC.getEnd(); pos++) {
                                        if(pos>=smallWinStart && pos<=smallWinEnd)
                                            bps.add(pos);
                                    }
                                }
                                regionVcs.add(lastVC);
                            }else if(curPosEndLong<smallWinStartLong){
                                remove=true;
                            }else{
                                breakFor=true;
                                break;
                            }

                        }
                        if(remove){
                            lastVCs.remove(i);
                        }
                    }
                    if(breakFor){
                        break;
                    }
                    String line = null;
                    while(true){
                        try {
                            if(samplesReader.get(i).getFilePointer()>endOffset.get(i)){
                                break;
                            }
                            line=samplesReader.get(i).readLine();
                            if(line==null){
                                readerDone[i]=true;
                                break;
                            }
                        } catch (IOException e) {
                            System.out.println("IO error!\tcurrent sample:\t"+samplesTag+"\t"+cycleIter);
                            System.out.println("current window:\t"+smallWinStartLong);
                            if(lastVCs.containsKey(i)) {
                                System.out.println("last VC:");
                                for(VariantContext lastVC:lastVCs.get(i)){
                                    System.out.println(lastVC);
                                }
                            }
                            e.printStackTrace();
                        }

                        VariantContext combineVC=samplesCodec.get(i).decode(line);
                        long curPosStartLong=dBC.accumulateLength.get(dBC.chrIndex.get(combineVC.getContig()))+combineVC.getStart();
                        long curPosEndLong=dBC.accumulateLength.get(dBC.chrIndex.get(combineVC.getContig()))+combineVC.getEnd();
                        if(curPosStartLong<=smallWinEndLong && curPosEndLong>=smallWinStartLong){
                            if(combineVC.getNAlleles()>2) {
                                for (int pos = combineVC.getStart(); pos <= combineVC.getEnd(); pos++) {
                                    if(pos>=smallWinStart && pos<=smallWinEnd)
                                        bps.add(pos);
                                }
                            }
                            regionVcs.add(combineVC);
                            if(curPosEndLong>smallWinEndLong){
                                if(lastVCs.containsKey(i)){
                                    lastVCs.get(i).add(combineVC);
                                }else{
                                    ArrayList<VariantContext> tmpVCs=new ArrayList<>();
                                    tmpVCs.add(combineVC);
                                    lastVCs.put(i,tmpVCs);
                                }
                                break;
                            }
                        }else if(curPosStartLong>smallWinEndLong){
                            if(lastVCs.containsKey(i)){
                                lastVCs.get(i).add(combineVC);
                            }else{
                                ArrayList<VariantContext> tmpVCs=new ArrayList<>();
                                tmpVCs.add(combineVC);
                                lastVCs.put(i,tmpVCs);
                            }
                            break;
                        }else{
                            continue;
                        }
                    }
                }
                regionVcs.sort(comparator4);
//                if(logIter%5==0) {
                    System.out.println(formatter.format(new Date()) + "\tcurrent reduce key:\t" + chr + "\t" + chrInx + "\t" + smallWinStartLong + "\t" + smallWinEndLong);
                    System.out.println(formatter.format(new Date()) + "\tbps size in region:\t" + bps.size());
//                }
                logIter++;
                for(int iter:bps){
                    VariantContext variantContext = engine.variantCallingForSpark(regionVcs.iterator(),
                            parser.createGenomeLocation(chr, iter), genomeShare.getChromosomeInfo(chr));
                    if (variantContext == null){
                        continue;
                    }
                    CommonInfo info = variantContext.getCommonInfo();
                    HashMap<String, Object> maps = new HashMap<>(info.getAttributes());
                    maps.remove("SM");
                    info.setAttributes(maps);
                    String value=vcfEncoder.encode(variantContext)+"\n";
                    long vcStart=dBC.accumulateLength.get(dBC.chrIndex.get(variantContext.getContig()))+variantContext.getStart();
                    boolean writeFlag=false;
                    if(index==0){
                        if(vcStart<bpPartition.get(index)){
                            writeFlag=true;
                        }
                    }else{
                        if(vcStart<bpPartition.get(index) && vcStart>=bpPartition.get(index-1)){
                            writeFlag=true;
                        }
                    }
                    if(writeFlag) {
                        out.write(value.getBytes());
                        totalVariantsNum.add(1);
                        vcNum++;
                    }
                    maps.clear();
                    while(regionVcs.size()>0){
                        VariantContext firstVc=regionVcs.getFirst();
                        if(firstVc.getEnd()<=iter){
                            regionVcs.removeFirst();
                        }else{
                            break;
                        }
                    }
                }
                smallWinStartLong=smallWinEndLong+1;
                smallWinEndLong=smallWinStartLong+dBC.options.getWindowsSize()-1;
                if(smallWinEndLong>regionEnd){
                    smallWinEndLong=regionEnd;
                }
                smallWinStart=smallWinEnd+1;
                smallWinEnd=smallWinStart+dBC.options.getWindowsSize()-1;
                if(smallWinEnd>contigLength){
                    smallWinEnd=contigLength;
                }
                if(smallWinEndLong-smallWinStartLong>contigLength || smallWinEnd-smallWinStart>contigLength){
                    log.error("code error");
                    System.exit(1);
                }
                regionVcs.clear();
                bps.clear();
            }
            if(dbsnps!=null)
                dbsnps.clear();
        }
        out.close();
        ArrayList<String> totalNum=new ArrayList<>();
        totalNum.add("done"+index);
        for(int i=0;i<samplesReader.size();i++){
            samplesReader.get(i).close();
        }
        return totalNum.iterator();
    }
}
