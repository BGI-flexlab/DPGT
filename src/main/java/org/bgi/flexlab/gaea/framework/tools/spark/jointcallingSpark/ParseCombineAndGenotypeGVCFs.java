package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
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
    private final MultipleVCFHeaderForJointCalling headers = new MultipleVCFHeaderForJointCalling();
    private final GenomeLongRegion region;
    private final SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
    private final HashMap<String, String> confMap;
    private final String bpFile;
    private final DriverBC dBC;
    private final ArrayList<GenomeLocation> regions;
    private final int cycleIter;
    private final ArrayList<Long> bpPartition=new ArrayList<>();
    private final LongAccumulator totalVariantsNum;
    public ParseCombineAndGenotypeGVCFs(GenomeLongRegion region, ArrayList<GenomeLocation> regions, String outputBP, HashMap<String, String> confMap,
                                        Broadcast<DriverBC> dBC, int cycleIter, ArrayList<Long> bpPartition,
                                        LongAccumulator totalVariantsNum) {
        this.dBC=dBC.value();
        this.cycleIter=cycleIter;
        this.confMap=confMap;
        this.region=region;
        this.regions=regions;
        bpFile=outputBP;
        this.totalVariantsNum=totalVariantsNum;
        this.bpPartition.addAll(bpPartition);
    }
    static class VcComp implements Comparator<VariantContext>,Serializable{

        @Override public int compare(VariantContext o1, VariantContext o2) {
            return o1.getStart()-o2.getStart();
        }
    }
    private final Comparator<VariantContext> comparator4 = new VcComp();

    @Override public Iterator<String> call(Integer index,Iterator<String> stringIterator) throws IOException {
        Logger log= LoggerFactory.getLogger(ParseCombineAndGenotypeGVCFs.class);
        Configuration hadoop_conf=new Configuration();
        for(String k:confMap.keySet()){
            hadoop_conf.set(k,confMap.get(k));
        }
        System.out.println(formatter.format(new Date())+"\tbefore readHeader");
        SeekableStream in = new SeekableFileStream(new File(dBC.outputDir+"/vcfheader"));
        VCFHeader header = VCFHeaderReader.readHeaderFrom(in);
        System.out.println(formatter.format(new Date())+"\tafter readHeader");
        in.close();
        if(header == null)
            throw new RuntimeException("header is null !!!");

        GenomeLocationParser parser = new GenomeLocationParser(header.getSequenceDictionary());
        hadoop_conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, dBC.options.getOutDir() + "/vcfheader");
        String MERGER_HEADER_INFO = "vcfheaderinfo";
        hadoop_conf.set(MERGER_HEADER_INFO, dBC.options.getOutDir()+"/"+ MERGER_HEADER_INFO);
        SeekableStream in2=new SeekableFileStream(new File(dBC.outputDir+"/vcfheader"));
        headers.setMergeHeader(VCFHeaderReader.readHeaderFrom(in2));
        headers.setCurrentIndex(dBC.sampleIndex.size());
        HashMap<Integer, Set<String>> tmpNameHeaders= new HashMap<>();
        for(Map.Entry<String,Integer> kv:dBC.sampleIndex.entrySet()){
            TreeSet<String> samples=new TreeSet<>();
            samples.add(kv.getKey());
            tmpNameHeaders.put(kv.getValue(),samples);
        }
        headers.setNameHeaders(tmpNameHeaders);
        int totalSampleSize=dBC.sampleIndex.size();
        if(totalSampleSize==0){
            log.error("sample size is 0");
            System.exit(1);
        }
        int realSmallWindowSize=totalSampleSize>100000?10000000/totalSampleSize:dBC.options.getWindowsSize();
        if(realSmallWindowSize==0){
            realSmallWindowSize=10;
        }
        String sampleStr = confMap.get(DriverBC.INPUT_ORDER);

        dBC.pathSample=null;
        System.out.println(formatter.format(new Date())+"\tbefore engine init");
        dBC.options.setS(dBC.options.STANDARD_CONFIDENCE_FOR_CALLING);
        dBC.options.sets(dBC.options.STANDARD_CONFIDENCE_FOR_EMITTING);
        JointCallingEngine engine = new JointCallingEngine(dBC.options, parser, header, headers, dBC.multiMapSampleNames, hadoop_conf);
        System.out.println(formatter.format(new Date())+"\tafter engine init");
        ReferenceShare genomeShare = new ReferenceShare();
        genomeShare.loadChromosomeList(dBC.options.getReference());
        DbsnpShare dbsnpShare = new DbsnpShare(dBC.options.getDBSnp(), dBC.options.getReference());
        dbsnpShare.loadChromosomeList(dBC.options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
        VCFLocalLoader loader = new VCFLocalLoader(dBC.options.getDBSnp());
        VariantRegionFilter filter = new VariantRegionFilter();
        VCFHeader header2 = engine.getVCFHeader();
        if(header2 == null)
            throw new RuntimeException("header is null !!!");
        String mergedHeaderFile=dBC.options.getOutDir()+"/merged.vcfheader";
        File hFile=new File(mergedHeaderFile);
        if(!hFile.exists()) {
            VCFLocalWriter mergedVCFHeaderWriter = new VCFLocalWriter(mergedHeaderFile, false, false);
            mergedVCFHeaderWriter.writeHeader(header2);
            mergedVCFHeaderWriter.close();
        }
        BufferedReader bp_reader=new BufferedReader(new FileReader(bpFile));
        String winLine = bp_reader.readLine();
        System.out.println("reduce setup done");

        System.out.println(formatter.format(new Date())+"\treduce start");
        ArrayList<BlockCompressedInputStream> samplesReader=new ArrayList<>();
        ArrayList<VCFCodec> samplesCodec=new ArrayList<>();
        ArrayList<Integer> samplesTag=new ArrayList<>();
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
            String targetFileName=dBC.options.getOutDir()+"/combine."+cycleIter+"/combineFile."+mapSMtagInt;
            System.out.println("processed file:\t"+targetFileName);

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
        dBC.sampleIndex=null;
        int multiMapSampleSize=dBC.multiMapSampleNames.size();
        dBC.multiMapSampleNames=null;
        if(samplesReader.size()!=multiMapSampleSize || samplesCodec.size()!=multiMapSampleSize){
            log.error("code error2\t"+samplesReader.size()+"\t"+multiMapSampleSize);
            System.exit(1);
        }
        String outFile=dBC.options.getOutDir()+"/genotype."+cycleIter+"/"+region.getStart()+"_"+region.getEnd()+"."+index+".vcf.gz";
        Path outPath=new Path(outFile);
        FileSystem outFs=outPath.getFileSystem(hadoop_conf);
        if(outFs.exists(outPath)){
            outFs.delete(outPath,true);
        }
        BlockCompressedOutputStream out=new BlockCompressedOutputStream(outFile);
        if (hFile.exists()) {
            BufferedReader hfReader = new BufferedReader(new FileReader(mergedHeaderFile));
            String hfLine;
            while ((hfLine = hfReader.readLine()) != null) {
                out.write((hfLine + "\n").getBytes());
            }
            hfReader.close();
        } else {
            log.error("no header file");
            System.exit(1);
        }

        long rangeInEachOutFile=(region.getEnd()-region.getStart())/dBC.options.getReducerNumber();
        log.info("range in each out file:\t"+rangeInEachOutFile);
        log.info("current index:\t"+index);

        boolean[] readerDone=new boolean[multiMapSampleSize];
        for(int i=0;i<multiMapSampleSize;i++){
            readerDone[i]=false;
        }
        VCFEncoder vcfEncoder = new VCFEncoder(header, true, true);
        for(GenomeLocation curRegion:regions) {
            int start = curRegion.getStart();
            String chr = curRegion.getContig();
            int chrInx = dBC.chrIndex.get(chr);

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
//
            //为了避免每次从头读取winBp，保留bp_reader，使其只读一遍
            System.out.println(formatter.format(new Date()) + "\tbp window before get bps:\t" + winLine);
            LinkedList<VariantContext> regionVcs=new LinkedList<>();
            long smallWinStartLong=regionStart;
            long smallWinEndLong=smallWinStartLong+realSmallWindowSize-1;
            if(smallWinEndLong>regionEnd){
                smallWinEndLong=regionEnd;
            }
            int smallWinStart=start;
            int smallWinEnd=smallWinStart+realSmallWindowSize-1;
            if(smallWinEnd>contigLength){
                smallWinEnd=contigLength;
            }
            if(smallWinEndLong-smallWinStartLong>contigLength || smallWinEnd-smallWinStart>contigLength){
                log.error("code error");
                System.exit(1);
            }
            HashMap<Integer,ArrayList<VariantContext>> lastVCs=new HashMap<>();
            while(smallWinStartLong<=regionEnd){
                //一个个小窗口开始处理
                //获取小窗口内的所有变异
                if(smallWinStartLong>bpPartition.get(index)){
                    break;
                }
                if(index>0){
                    if(smallWinEndLong<bpPartition.get(index-1)) {
                        smallWinStartLong=smallWinEndLong+1;
                        smallWinEndLong=smallWinStartLong+realSmallWindowSize-1;
                        if(smallWinEndLong>regionEnd){
                            smallWinEndLong=regionEnd;
                        }
                        smallWinStart=smallWinEnd+1;
                        smallWinEnd=smallWinStart+realSmallWindowSize-1;
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
                int minWindowSize=dBC.options.getWindowsSize();
                if(realSmallWindowSize<100){
                    minWindowSize=100;
                }
                long startPosition = dbsnpShare.getStartPosition(chr,smallWinStart / minWindowSize,smallWinEnd/minWindowSize,minWindowSize);
                if (startPosition >= 0) dbsnps = filter.loadFilter(loader, chr, startPosition, smallWinEnd);
                engine.init(dbsnps);
                System.out.println("dbsnp size:\t"+dbsnps.size());
                Set<Integer> bps=new TreeSet<>();
                for(int i=0;i<multiMapSampleSize;i++){
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
                        }
                    }
                }
                regionVcs.sort(comparator4);
                System.out.println(formatter.format(new Date()) + "\tcurrent reduce key:\t" + chr + "\t" + chrInx + "\t" + smallWinStartLong + "\t" + smallWinEndLong);
                System.out.println(formatter.format(new Date()) + "\tbps size in region:\t" + bps.size());
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
                    String value= vcfEncoder.encode(variantContext)+"\n";
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
                smallWinEndLong=smallWinStartLong+realSmallWindowSize-1;
                if(smallWinEndLong>regionEnd){
                    smallWinEndLong=regionEnd;
                }
                smallWinStart=smallWinEnd+1;
                smallWinEnd=smallWinStart+realSmallWindowSize-1;
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
            dbsnps.clear();
        }
        out.close();
        ArrayList<String> totalNum=new ArrayList<>();
        totalNum.add("done"+index);
        for (BlockCompressedInputStream blockCompressedInputStream : samplesReader) {
            blockCompressedInputStream.close();
        }
        return totalNum.iterator();
    }
}
