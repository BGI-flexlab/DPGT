package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.broadcast.Broadcast;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.variant.VariantContextMerger;
import org.bgi.flexlab.gaea.tools.jointcalling.VariantAnnotatorEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.ReferenceConfidenceVariantContextMerger;
import org.bgi.flexlab.gaea.util.Utils;
import org.seqdoop.hadoop_bam.LazyVCFGenotypesContext;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;


public class CombineVariants implements VoidFunction<Iterator<String>> {
    private final GenomeLongRegion gloc;
    private final String bpPath;
    private final String vVcfPath;
//    public Set<VCFHeaderLine> gvcfHeaderMetaInfo;
//    public VCFHeaderVersion version = null;
private JointCallingSparkOptions options = new JointCallingSparkOptions();
    private final static String INPUT_LIST = "input.gvcf.list";
    private HashMap<String, String> confMap=new HashMap<>();
    private final String outputDir;
    private final Logger logger;
    private VCFEncoder vcfEncoder=null;
    private final DriverBC dBC;
    private ArrayList<GenomeLocation> regions=new ArrayList<>();
    private final int cycleIter;
    private final ArrayList<Long> bpPartition=new ArrayList<>();
//    public static VCFHeader mergedHeader=null;
    public CombineVariants(GenomeLongRegion region, ArrayList<GenomeLocation> regions,String outputBP, HashMap<String, String> confMap,
                           Broadcast<DriverBC> dBC,int cycleIter,ArrayList<Long> bpPartition) {
        gloc=region;
        this.regions=regions;
        bpPath=outputBP;
        this.cycleIter=cycleIter;
        this.dBC=dBC.value();
        this.options=dBC.value().options;
        this.confMap=confMap;
        this.outputDir=dBC.value().outputDir;
        this.confMap.put(INPUT_LIST, options.getInputList());
        logger = LoggerFactory.getLogger(CombineVariants.class);
        this.vVcfPath=dBC.value().vVcfPath;
        this.bpPartition.addAll(bpPartition);

    }
    class VcComp implements Comparator<VariantContext>,Serializable{

        @Override public int compare(VariantContext o1, VariantContext o2) {
            return o1.getStart()-o2.getStart();
        }
    }
    private final Comparator<VariantContext> comparator4 = new VcComp();
    @Override public void call(Iterator<String> StringIterator) throws Exception {
        //与map的setup功能类似
        Configuration conf=new Configuration();
        SeekableStream in2 = new SeekableFileStream(new File(outputDir+"/vcfheader"));
        VCFHeader mergedHeader = VCFHeaderReader.readHeaderFrom(in2);
        if(mergedHeader==null){
            logger.error("header is null, "+in2.getSource());
            System.exit(1);
        }

        Configuration hadoop_conf=new Configuration();
        for(String k:confMap.keySet()){
            hadoop_conf.set(k,confMap.get(k));
        }
        ArrayList<VariantContext> vcList=new ArrayList<>();
        Set<String> mapSamplePath=new TreeSet<>();
        Set<VCFHeaderLine> gvcfHeaderMetaInfo=null;
        HashMap<String, Integer> sampleIndex=new HashMap();
        String sampleIndexFile=options.getOutDir()+"/vcfheaderinfo";
        File sampleIndexFilePath=new File(sampleIndexFile);
        BufferedReader indexReader = new BufferedReader(new FileReader(sampleIndexFilePath));
        String indexLine;
        Logger logger= LoggerFactory.getLogger(CombineVariants.class);
        Map<String,String> pathSample=new HashMap();
        ArrayList<String> sampleNames=new ArrayList<String>();
        Integer mapSMtagInt=0;
        int totalSampleSize=0;
        VCFHeaderVersion version = null;
//        HashMap<Integer,Set<String>> mapIndexSample=new HashMap<>();
//        ArrayList<VCFCodec> codec = new ArrayList<VCFCodec>();
        ArrayList<String> mapGvcfList=new ArrayList<String>();
        while(StringIterator.hasNext()) {
            String gvcfPath = StringIterator.next();
            mapGvcfList.add(gvcfPath);
        }

        for(String gvcfPath:mapGvcfList) {
//            Path path2=new Path(gvcfPath);
            String[] eles=gvcfPath.split("/");
            mapSamplePath.add(eles[eles.length-1]);
        }
        while((indexLine=indexReader.readLine())!=null) {
            totalSampleSize++;
            String[] eles=indexLine.split("\t");
            if(eles.length!=3) {
                logger.error("vcfheaderinfo file format error");
            }
            String name;
            if(eles[2].endsWith(",")) {
                name=eles[2].substring(0,eles[2].length()-1);
            }else {
                name=eles[2];
            }
            if(mapSamplePath.contains(eles[0])) {
                sampleIndex.put(name, Integer.parseInt(eles[1]));
                pathSample.put(eles[0], name);
            }
        }
        ArrayList<String> samplesInMap=new ArrayList<>();
        for(String samplePath:mapSamplePath){
            samplesInMap.add(pathSample.get(samplePath));
        }
        indexReader.close();
        HashMap<String,Integer> chrIndexs = null;
        VCFHeader virtualHeader = null;
        VCFHeader mapMergedHeader=null;
        String mapSMtag=null;
//        HashMap<Integer, String> contigs = new HashMap<Integer, String>();
        Map<String, Integer> contigDict = new HashMap<String, Integer>();
        GenomeLocationParser parser = null;
        ReferenceShare genomeShare = null;
//        DbsnpShare dbsnpShare = null;
        LazyVCFGenotypesContext.HeaderDataCache vcfHeaderDataCache =
                new LazyVCFGenotypesContext.HeaderDataCache();
        SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
        VariantAnnotatorEngine annotationEngine=null;
//        List<String> annotationsToUse = new ArrayList<>();

//        List<String> annotationGroupsToUse = new ArrayList<>(
//                Arrays.asList(new String[] { StandardAnnotation.class.getSimpleName() }));

//        ArrayList<BufferedReader> gvcfBufReader=new ArrayList<BufferedReader>();
        ArrayList<VariantContext> curSamplesVC=new ArrayList<VariantContext>();

        Integer sISize=sampleIndex.size();
        int ii=0;
        ArrayList<TabixFeatureReader> samplesReader=new ArrayList<>();
        SeekableStream in=new SeekableFileStream(new File(options.getOutDir()+"/virtual.vcf"));
        virtualHeader = VCFHeaderReader.readHeaderFrom(in);//从header文件中读取header
        in.close();
        VCFHeader curSamplesMergedHeader=new VCFHeader(virtualHeader.getMetaDataInInputOrder(),samplesInMap);
        if(vcfEncoder==null) {
            vcfEncoder = new VCFEncoder(curSamplesMergedHeader, true, false);
        }
        mapSMtagInt=totalSampleSize;
        int mapGvcfListIndex=0;
        gvcfHeaderMetaInfo=virtualHeader.getMetaDataInInputOrder();
//        Set<String> chrs=new TreeSet<>();
//        for(GenomeLocation l:regions){
//            chrs.add(l.getContig());
//        }
        for(String gvcfPath:mapGvcfList){
            String[] eles=gvcfPath.split("/");
            mapSamplePath.add(eles[eles.length-1]);
            Path path2=new Path(gvcfPath);
            if(!pathSample.containsKey(eles[eles.length-1])) {
                logger.error("no such path in vcfHeaderInfo");
            }
            String sampleName=pathSample.get(eles[eles.length-1]);
            sampleNames.add(sampleName);
            if(mapGvcfListIndex==0) {
                mapSMtagInt += sampleIndex.get(sampleName);
            }
            mapGvcfListIndex++;
            if(version==null) {
                for (final VCFHeaderLine line : gvcfHeaderMetaInfo) {
                    if (VCFHeaderVersion.isFormatString(line.getKey())) {
                        version = VCFHeaderVersion.toHeaderVersion(line.getValue());
                        break;
                    }
                }
            }
            VCFCodec tmp_codec=new VCFCodec();
            Set<String> curSample=new HashSet<>();
            curSample.add(sampleName);
//            mapIndexSample.put(ii,curSample);
            tmp_codec.setVCFHeader(new VCFHeader(gvcfHeaderMetaInfo,curSample), version);
//            codec.add(tmp_codec);
            logger.warn(ii+"\tcurrent index");
            if(gvcfPath.startsWith("file://")) {
                gvcfPath=gvcfPath.substring(7);
            }
//            VCFCodec query_codec=new VCFCodec();
            TabixFeatureReader sampleReader=new TabixFeatureReader(gvcfPath,tmp_codec);
//            Iterator<VariantContext> it=sampleReader.query(gloc.getContig(),gloc.getStart(),gloc.getEnd());
            samplesReader.add(sampleReader);
//            if(gvcfPath.endsWith(".gz")) {
//                if(gvcfPath.startsWith("/user")) {
//                    reader=new BufferedReader(new InputStreamReader(new GZIPInputStream(path2.getFileSystem(hadoop_conf).open(path2))));
//                    reader.skip(1);
//                }else {
//                    reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gvcfPath))));
//                }
//            }else {
//                reader = new BufferedReader(new InputStreamReader(new FileInputStream(gvcfPath)));
//            }
//            gvcfBufReader.add(reader);
            ii++;

        }
        gvcfHeaderMetaInfo=virtualHeader.getMetaDataInInputOrder();
        mapSamplePath.clear();


        logger.warn("sampleNames Size:\t"+sampleNames.size());
        logger.warn("mapGvcfList Size:\t"+mapGvcfList.size());
        pathSample.clear();
        sampleIndex.clear();

        long rangeInEachOutFile=(gloc.getEnd()-gloc.getStart())/options.getReducerNumber();
        int minRange=100000;
        if(rangeInEachOutFile<minRange){
            rangeInEachOutFile=minRange;
        }
        String oneOutFile=options.getOutDir()+"/combine."+cycleIter+"/combineFile."+mapSMtagInt;
        File outOutFilePath=new File(oneOutFile);
        if(outOutFilePath.exists()){
            outOutFilePath.delete();
        }
//        BufferedWriter outWriter=new BufferedWriter(new FileWriter(outFile));
//        long curOffset=0;
//        TreeMap<Integer,Long> outputIndexOffset=new TreeMap<>();
        BlockCompressedOutputStream outWriter=new BlockCompressedOutputStream(oneOutFile);
//        ArrayList<BufferedWriter> outWriterList=new ArrayList<>();
        for(int i=0;i<options.getReducerNumber();i++){
            String outFile=options.getOutDir()+"/combine."+cycleIter+"/combineFile."+mapSMtagInt+"."+i;
//            Path outPath=new Path(outFile);
//            FileSystem outFs=outPath.getFileSystem(hadoop_conf);
//            if(outFs.exists(outPath)){
//                outFs.delete(outPath,true);
//            }
//            BufferedWriter curRegionTheseSamples=new BufferedWriter(new OutputStreamWriter(outFs.create(outPath)));
            File outFilePath=new File(outFile);
            if(outFilePath.exists()){
                outFilePath.delete();
            }
//            BufferedWriter curRegionTheseSamples=new BufferedWriter(new FileWriter(outFile));
//            outWriterList.add(curRegionTheseSamples);
        }


//        System.out.println("mapSMtag:\t"+mapSMtagInt);
        mapMergedHeader=new VCFHeader(gvcfHeaderMetaInfo,sampleNames);
        //mapMergedHeader=getVCFHeaderFromInput(mapMultiSamplesHeaderSet);
        vcfHeaderDataCache.setHeader(mapMergedHeader);

        mapSMtag= Utils.join("_",mapMergedHeader.getSampleNamesInOrder());
        if(mapMergedHeader.getSampleNamesInOrder().size()==1) {
            mapSMtag=mapSMtag+"_";
        }
        if(virtualHeader == null)
            throw new RuntimeException("header is null !!!");
        logger.warn("after get merged header");
//        contigs = new HashMap<>();
        chrIndexs = new HashMap<>();
        for (VCFContigHeaderLine line : virtualHeader.getContigLines()) {
            chrIndexs.put(line.getID(), line.getContigIndex());
        }
//        for (VCFContigHeaderLine line : virtualHeader.getContigLines()) {
//            contigs.put(line.getContigIndex(), line.getID());
//        }

        parser = new GenomeLocationParser(virtualHeader.getSequenceDictionary());

        logger.warn("after parser");

        genomeShare = new ReferenceShare();
        genomeShare.loadChromosomeList(options.getReference());
//        dbsnpShare = new DbsnpShare(options.getDBSnp(), options.getReference());
//        dbsnpShare.loadChromosomeList(options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
        for (final VCFContigHeaderLine contig : virtualHeader.getContigLines())
            contigDict.put(contig.getID(), contig.getContigIndex());
        logger.warn("after contigDict Map done");
//        annotationEngine = new VariantAnnotatorEngine(annotationGroupsToUse, annotationsToUse,
//                Collections.<String>emptyList());
//        annotationEngine.initializeDBs(options.getDBSnp() != null);
        logger.warn("setup done");

        //与map的map功能类似
        String winFilePath=hadoop_conf.get(dBC.Window_File);
        if(winFilePath.startsWith("file://")) {
            winFilePath=winFilePath.substring(7);
        }
        String win_line;
//        gvcfBufReader.clear();
        Path firstPath=new Path(mapGvcfList.get(0));
        FileSystem gvcf_fs=firstPath.getFileSystem(hadoop_conf);
//        for(int i=0;i!=mapGvcfList.size();i++) {
//            String gvcfFile=mapGvcfList.get(i);
//            if(gvcfFile.startsWith("file://")) {
//                gvcfFile=gvcfFile.substring(7);
//            }
//            BufferedReader reader=null;
//            if(gvcfFile.endsWith(".gz")) {
//                if(gvcfFile.startsWith("/user")) {
//                    reader=new BufferedReader(new InputStreamReader(new GZIPInputStream(gvcf_fs.open(new Path(gvcfFile)))));
//                }else {
//                    reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gvcfFile))));
//                }
//            }else {
//                reader = new BufferedReader(new InputStreamReader(new FileInputStream(gvcfFile)));
//            }
//            gvcfBufReader.add(reader);
//        }
//        String[] tmpLines=new String[mapGvcfList.size()];
//        for(int i=0;i!=mapGvcfList.size();i++) {
//            String tmpline=gvcfBufReader.get(i).readLine();
//            tmpLines[i]=tmpline;
//        }
        //method 3, merge code from MapperHuge
        BufferedReader win_reader=new BufferedReader(new FileReader(winFilePath));
        String last_contig="";
        ChromosomeInformationShare ref=null;
        Path pBPpath=new Path(bpPath);
        BufferedReader bp_reader=new BufferedReader(new FileReader(bpPath));
        String bpRegion=bp_reader.readLine();
        if(bpRegion==null){
            return;
        }
        String[] bpeles=bpRegion.split("\t");
        Long bpstart=Long.parseLong(bpeles[0]);
        Long bpend=Long.parseLong(bpeles[1]);
        ArrayList<Integer> bpStarts=new ArrayList<>();
        ArrayList<Integer> bpEnds=new ArrayList<>();
        boolean bpFileEnd=false;
        boolean logOut=false;
        Integer lastLogPoint=0;
        Integer last_end=0;
        Integer last_pos=0;
        Long last_bp_start=-1L;
        Long last_bp_end=-1L;
//        Set<Integer> last_end_breakpoints=new TreeSet<Integer>();
        ArrayList<Iterator<VariantContext>> sampleVCsIter=new ArrayList<>();
//        for(int i=0;i!=mapGvcfList.size();i++) {
//            VCFFileReader fileReader = dBC.sampleReader.get(mapGvcfList.get(i));
//            Iterator<VariantContext> vcsAtRegion = fileReader.query(gloc.getContig(), gloc.getStart(), gloc.getEnd());
//            sampleVCsIter.add(vcsAtRegion);
//        }

        ArrayList<Iterator<VariantContext>> samplesIt=new ArrayList<>();
        ArrayList<ArrayList<VariantContext>> samplesVcList=new ArrayList<>();
        int[] VcIndex=new int[mapGvcfList.size()];
        TreeMap<Integer,Long> startOffset=new TreeMap<>();
        TreeMap<Integer,Long> endOffset=new TreeMap<>();
        logger.info("region size:\t"+regions.size());
        for(GenomeLocation curRegion:regions){
            logger.info(curRegion.toString());
        }
        for(GenomeLocation curRegion:regions) {
            String contig=curRegion.getContig();
            // start end 为小窗口
            Integer start=curRegion.getStart();
            Integer end=curRegion.getEnd()>start+dBC.options.getWindowsSize()-1?start+dBC.options.getWindowsSize()-1:curRegion.getEnd();
            Long smallWinStart=dBC.accumulateLength.get(dBC.chrIndex.get(curRegion.getContig()))+start;
            Long smallWinEnd=dBC.accumulateLength.get(dBC.chrIndex.get(curRegion.getContig()))+end;
            Long regionEnd=dBC.accumulateLength.get(dBC.chrIndex.get(curRegion.getContig()))+curRegion.getEnd();
            for (int i = 0; i != mapGvcfList.size(); i++) {
                Iterator<VariantContext> curIt = samplesReader.get(i).query(curRegion.getContig(), curRegion.getStart(), curRegion.getEnd());
//            ArrayList<VariantContext> tmpVCs=new ArrayList<>();
//            while(curIt.hasNext()){
//                tmpVCs.add(curIt.next());
//            }
//            samplesVcList.add(tmpVCs);
                logger.info("process samples in this partition:\t" + mapGvcfList.get(i));
                samplesIt.add(curIt);
//            VcIndex[i]=0;
//            samplesReader.get(i).close();
            }
            HashMap<Integer, VariantContext> lastVC = new HashMap<>();
            int testInt = 0;
            while (smallWinEnd <= regionEnd) {
                Integer curLogPoint = (int) (smallWinStart / 1000000);
                if (!curLogPoint.equals(lastLogPoint)) {
                    logOut = true;
                } else {
                    logOut = false;
                }
                lastLogPoint = curLogPoint;
                if (logOut)
                    System.out.println(formatter.format(new Date()) + "\tprocess region\t" + contig + "\t" + smallWinStart + "\t" + smallWinEnd);
                String chr = contig;
                int startPartitionIndex=0;
                for(int partitionIndex=0;partitionIndex<bpPartition.size();partitionIndex++){
                    if(smallWinStart<bpPartition.get(partitionIndex)){
                        startPartitionIndex=partitionIndex;
                        break;
                    }
                }
                Integer winChrID = chrIndexs.get(contig);
//                if (last_bp_start != -1) {
////                    bpStarts.add(last_bp_start);
////                    bpEnds.add(last_bp_end);
////                    last_bp_start = -1L;
////                    last_bp_end = -1L;
//                }
                //将目标处理区域内的bp信息提取出来
                while (!bpFileEnd) {
                    if (bpstart <= smallWinEnd && bpend >= smallWinStart) {
                        int intBpStart= (int) (bpstart-dBC.accumulateLength.get(winChrID));
                        int intBpEnd= (int) (bpend-dBC.accumulateLength.get(winChrID));

                        bpStarts.add(intBpStart);
                        bpEnds.add(intBpEnd);
                        if(bpend>smallWinEnd){
//                            last_bp_start = bpstart;
//                            last_bp_end = bpend;
                            break;
                        }
                        bpRegion = bp_reader.readLine();
                        if (bpRegion == null) {
                            bpFileEnd = true;
                            bp_reader.close();
                            break;
                        }
                        bpeles=bpRegion.split("\t");
                        bpstart=Long.parseLong(bpeles[0]);
                        bpend=Long.parseLong(bpeles[1]);
                    } else if (bpstart > smallWinEnd) {
//                        last_bp_start = bpstart;
//                        last_bp_end = bpend;
                        break;
                    } else {
                        bpRegion = bp_reader.readLine();
                        if (bpRegion == null) {
                            bpFileEnd = true;
                            bp_reader.close();
                            break;
                        }
                        bpeles=bpRegion.split("\t");
                        bpstart=Long.parseLong(bpeles[0]);
                        bpend=Long.parseLong(bpeles[1]);
                    }
                }
                LinkedList<VariantContext> region_vcs = new LinkedList<>();
                Set<Integer> end_breakpoints = new TreeSet<>();
                if (logOut) {
                    System.out.println(formatter.format(new Date()) + "\tbefore query");
                }
                for (int i = 0; i != mapGvcfList.size(); i++) {
                    String sampleName = sampleNames.get(i);
                    //非query方式的代码起始位置
                /*
                VariantContext lastVc=null;
                File lastVCFile=new File((new File(vVcfPath)).getParent()+"/"+sampleName+".mr2.lastVC");
                if(lastVCFile.exists()){
                    BufferedReader lastVCReader=new BufferedReader(new FileReader(lastVCFile));
                    String line1Sample,line2VC;
                    if((line1Sample=lastVCReader.readLine())!=null){
                        if(!line1Sample.equals(sampleName)){
                            logger.error("code error");
                            System.exit(1);
                        }
                    }
                    if((line2VC=lastVCReader.readLine())!=null){
                        lastVc=codec.get(i).decode(line2VC);
                    }
                    lastVCReader.close();
                }
                if(lastVc!=null){
                    if(lastVc.getChr().equals(chr) && lastVc.getStart()<=end && lastVc.getEnd()>=start) {
                        region_vcs.add(lastVc);
                        if(lastVc.getNAlleles()>2) {
                            for(int tmpPos=lastVc.getStart();tmpPos<=lastVc.getEnd();tmpPos++) {
                                end_breakpoints.add(tmpPos);
                            }
                        }else {
                            end_breakpoints.add(lastVc.getEnd());
                        }
                    }
                }
                */
//                    if (lastVC.containsKey(i)) {
//                        VariantContext lastVc = lastVC.get(i);
//                        if (lastVc.getChr().equals(chr) && lastVc.getStart() <= end && lastVc.getEnd() >= start) {
//                            region_vcs.add(lastVc);
//                            if (lastVc.getNAlleles() > 2) {
//                                for (int tmpPos = lastVc.getStart(); tmpPos <= lastVc.getEnd(); tmpPos++) {
//                                    end_breakpoints.add(tmpPos);
//                                }
//                            } else {
//                                end_breakpoints.add(lastVc.getEnd());
//                            }
//                        }
//                        lastVC.remove(i);
//                    }
                    //非query方式的代码结束位置
                }
                for (int i = 0; i != curSamplesVC.size(); i++) {
                    VariantContext tmpVC = curSamplesVC.get(i);
                    if (tmpVC.getChr().equals(chr) && tmpVC.getStart() <= end && tmpVC.getEnd() >= start) {
                        region_vcs.add(tmpVC);
                        if (tmpVC.getNAlleles() > 2) {
                            for (int tmpPos = tmpVC.getStart(); tmpPos <= tmpVC.getEnd(); tmpPos++) {
                                end_breakpoints.add(tmpPos);
                            }
                        } else {
                            end_breakpoints.add(tmpVC.getEnd());
                        }
                        if (tmpVC.getEnd() <= end) {
                            curSamplesVC.remove(i);
                            i--;
                        }
                    }
                }
                int bpStartsIndex=0;
                int bpEndsIndex=0;
                if(bpStarts.size()!=bpEnds.size()){
                    logger.error("code error");
                    System.exit(1);
                }
                for (int i = 0; i != mapGvcfList.size(); i++) {
                    String line = null;
                    String samplePath = mapGvcfList.get(i);

                    String sampleName = sampleNames.get(i);
//                Integer sampleIdx=dBC.sampleIndex.get(sampleName);
//                while((line=JointCallingSpark.sampleReadersForMR2.get(sampleIdx).readLine())!=null) {
//                Iterator<VariantContext> curIt=samplesReader.get(i).query(contig,start,end);
//                while(samplesVcList.get(i).size()>VcIndex[i]){
                    int index=0;
                    while (samplesIt.get(i).hasNext()) {
                        final VariantContext v = samplesIt.get(i).next();
//                    if(line.startsWith("#")) {
//                        continue;
//                    }
//                    final VariantContext v=samplesVcList.get(i).get(VcIndex[i]);
//                    VcIndex[i]++;

//                    final VariantContext v = codec.get(i).decode(line);
                        CommonInfo info = v.getCommonInfo();
                        if (!info.hasAttribute("SM")) info.putAttribute("SM", sampleName);
                        if (v.getChr().equals(chr)) {
                            if (v.getStart() <= end && v.getEnd() >= start) {
                                while(index<bpStarts.size()) {
                                    if(v.getStart()<=bpEnds.get(index) && v.getEnd()>=bpStarts.get(index)) {
                                        if (v.getNAlleles() > 2) {
                                            for (int tmpPos = v.getStart(); tmpPos <= v.getEnd(); tmpPos++) {
                                                end_breakpoints.add(tmpPos);
                                            }
                                        } else {
                                            end_breakpoints.add(v.getEnd());
                                        }
                                        region_vcs.add(v);
                                        break;
                                    }else if(v.getStart()>bpEnds.get(index)){
                                        index++;
                                    }else{
                                        break;
                                    }
                                }
                                if (v.getEnd() > end) {
                                    curSamplesVC.add(v);
                                    break;
                                }
                            } else if (v.getStart() > end) {
                                curSamplesVC.add(v);
//                                lastVC.put(i, v);
                                break;
                            } else if (v.getEnd() < start) {
                                continue;
                            } else {
                            }
                        } else if (contigDict.get(v.getChr()) > contigDict.get(chr)) {
                            curSamplesVC.add(v);
//                            lastVC.put(i, v);
                            break;
//                        writeToDisk(sampleName,v);
                        } else {
                            continue;
                        }
                    }
                }
                if (logOut) {
                    System.out.println(formatter.format(new Date()) + "\tafter query");
                }
                region_vcs.sort(comparator4);
                //do combineGVCFS

                if (!contig.equals(last_contig)) {
                    ref = genomeShare.getChromosomeInfo(contig);
                    last_pos = 0;
                } else {
                    last_pos = start;
                }
                ArrayList<Integer> realBreakpointStart = new ArrayList<Integer>();
                ArrayList<Integer> realBreakpointEnd = new ArrayList<Integer>();
                int bpIndex = 0;
                for (Integer pos : end_breakpoints) {
                    boolean whether_keep = false;
                    for (; bpIndex != bpStarts.size(); bpIndex++) {
                        if (last_pos <= bpEnds.get(bpIndex) && pos >= bpStarts.get(bpIndex)) {
                            whether_keep = true;
                            break;
                        } else if (last_pos > bpEnds.get(bpIndex)) {
                            continue;
                        } else if (pos < bpStarts.get(bpIndex)) {
                            break;
                        }
                    }
                    if (whether_keep) {
                        realBreakpointStart.add(last_pos);
                        realBreakpointEnd.add(pos);
                        if (pos > end) {
//                            last_end_breakpoints.add(pos);
                        }
                    }
                    last_pos = pos + 1;
                }

                if (region_vcs.size() > 0) {
                    if (logOut) {
                        System.out.println(formatter.format(new Date()) + "\tvc number in this region:\t" + region_vcs.size() + "\tcurSamplesVC:\t" + curSamplesVC.size() + "\tbpStarts:\t" + bpStarts.size() + "\tend_breakpoints:\t" + end_breakpoints.size() + "\trealBreakpoints:\t" + realBreakpointStart.size());
                    }
                }
                end_breakpoints.clear();
                bpStarts.clear();
                bpEnds.clear();
                if (!contig.equals(last_contig)) {
                    last_end = 0;
                }
                //last_pos=0;
                int tmpIter = 0;
                if (logOut) {
                    System.out.println(formatter.format(new Date()) + "\tbefore merge");
                }
                for (int i = 0; i != realBreakpointStart.size(); i++) {
                    Integer posStart = realBreakpointStart.get(i);
                    Integer pos = realBreakpointEnd.get(i);

                    if (pos > end) {
                        break;
                    }
                    if (pos < start) {
                        continue;
                    }
                    List<VariantContext> stoppedVCs = new ArrayList<>();
                    //System.out.println(region_vcs.size());
                    for (int j = 0; j < region_vcs.size(); j++) {
                        final VariantContext vc = region_vcs.get(j);
                        //the VC for the previous state will be stopped if its position is previous to the current position or it we've moved to a new contig
                        if (vc.getStart() <= pos) {

                            if (vc.getEnd() < posStart) {
                                region_vcs.remove(j);
                                j--;
                            } else {
                                stoppedVCs.add(vc);
                            }
                            // if it was ending anyways, then remove it from the future state
                            if (vc.getEnd() == pos) {
                                //region_vcws.samples.removeAll(vc.getSampleNames());
                                if (j >= 0) {
                                    region_vcs.remove(j);
                                    j--;
                                }
                            }
                        } else {
                            break;
                        }
                    }
                    //System.out.println(region_vcs.size());
                    if (!stoppedVCs.isEmpty()) {
                        //System.out.println(stoppedVCs.size());
//					for(VariantContext vc:stoppedVCs) {
//						System.out.println(vc);
//					}
                        final GenomeLocation gLoc = parser.createGenomeLocation(chr, last_end);
                        final VariantContext mergedVC;
                        if(posStart>=dBC.accumulateLength.get(gLoc.getContigIndex()+1)){
                            break;
                        }
                        final byte refBase = (byte) ref.getBase(gLoc.getStart());
                        int curStart = last_end + 1;
                        if (containsTrueAltAllele(stoppedVCs)) {
                            mergedVC = ReferenceConfidenceVariantContextMerger.mapMerge(stoppedVCs, parser.createGenomeLocation(chr, pos), (byte) ref.getBase(posStart - 1), false, false, annotationEngine);
                            tmpIter++;
//						if(mergedVC.getStart()==1002417){
//							System.out.println(mergedVC+"\n"+mergedVC.getAttribute("DP"));
//						}
                        } else {
                            mergedVC = referenceBlockMerge(stoppedVCs, gLoc, refBase, pos);
//						if(mergedVC.getStart()==1002417){
//							System.out.println(mergedVC+"\n"+mergedVC.getAttribute("DP"));
//						}
                        }
                        if (mergedVC == null) {
                            continue;
                        }

                        CommonInfo info = mergedVC.getCommonInfo();
                        if (!info.hasAttribute("SM")){
                            info.putAttribute("SM", mapSMtagInt);
                        }else{
                            logger.warn("SM tag exists"+info.getAttribute("SM")+"\t"+mergedVC);
                            info.putAttribute("SM", mapSMtagInt);
                        }
                        if(!info.hasAttribute("END")){
                            info.putAttribute("END",mergedVC.getEnd());
                        }
                        HashMap<String, Object> maps = new HashMap<>(info.getAttributes());
                        info.setAttributes(maps);
//                        GenomeLocation vcLoc = new GenomeLocation(mergedVC.getContig(), chrIndexs.get(mergedVC.getContig()), mergedVC.getStart(), mergedVC.getEnd());
                        //有返回值的写法
//                    vcList.add(mergedVC);
                        long vcStartLong=dBC.accumulateLength.get(chrIndexs.get(mergedVC.getContig()))+mergedVC.getStart();
                        long vcEndlong=dBC.accumulateLength.get(chrIndexs.get(mergedVC.getContig()))+mergedVC.getEnd();
                        int startIter=0;
                        int endIter=0;
                        for(int partitionIndex=startPartitionIndex;partitionIndex<bpPartition.size();partitionIndex++){
                            if(vcStartLong<bpPartition.get(partitionIndex)){
                                startIter=partitionIndex;
                            }
                            if(vcEndlong<bpPartition.get(partitionIndex)){
                                endIter=partitionIndex;
                                break;
                            }
                        }
//                        int startIter= (int) (vcStartLong<gloc.getStart()?0:(int)((vcStartLong-gloc.getStart())/rangeInEachOutFile));
//                        int endIter= vcEndlong>gloc.getEnd()?options.getReducerNumber()-1: (int) ((vcStartLong - gloc.getStart()) / rangeInEachOutFile);
                        if(startIter>=options.getReducerNumber()){
                            startIter=options.getReducerNumber()-1;
                        }
                        if(endIter>=options.getReducerNumber()){
                            endIter=options.getReducerNumber()-1;
                        }
                        if(startIter<endIter){
                            logger.error("code error");
                            System.exit(1);
                        }


                        for(int iter=startIter;iter<=endIter;iter++) {
                            if(!startOffset.containsKey(iter)){
                                startOffset.put(iter,outWriter.getFilePointer());
                            }
                            endOffset.put(iter,outWriter.getFilePointer());
//                            String vcLine=vcfEncoder.encode(mergedVC)+"\n";
//                            logger.info("vc start:\t"+mergedVC.getStart()+", will write to "+mapSMtagInt+"."+iter+"\t"+rangeInEachOutFile+"\t"+gloc.getStart());

//                            outWriterList.get(iter).write(vcLine+"\n");
//                            System.out.println(outWriter.getFilePointer());

                        }
                        String vcLine=vcfEncoder.encode(mergedVC)+"\n";
                        outWriter.write(vcLine.getBytes());
                        stoppedVCs.clear();
                    }
                    last_end = pos;
                }
                if (logOut) {
                    System.out.println(formatter.format(new Date()) + "\tafter merge");
                }
                last_contig = contig;

                end_breakpoints.clear();
                start = end + 1;
                end = curRegion.getEnd() > end + dBC.options.getWindowsSize() ? end + dBC.options.getWindowsSize() : curRegion.getEnd();
                smallWinStart=dBC.accumulateLength.get(dBC.chrIndex.get(curRegion.getContig()))+start;
                smallWinEnd=dBC.accumulateLength.get(dBC.chrIndex.get(curRegion.getContig()))+end;
                if (start > end) {
                    break;
                }
            }
            samplesIt.clear();
        }
        win_reader.close();
//        for(int i=0;i!=outWriterList.size();i++){
//            outWriterList.get(i).close();
//        }
        outWriter.close();
        File idxFile=new File(oneOutFile+".idx");
        if(idxFile.exists()){
            idxFile.delete();
        }
        BufferedWriter idxWriter=new BufferedWriter(new FileWriter(oneOutFile+".idx"));
        for(Map.Entry<Integer,Long> kv:startOffset.entrySet()){
            if(!endOffset.containsKey(kv.getKey())){
                logger.error("code error");
                System.exit(1);
            }
            idxWriter.write(kv.getKey()+"\t"+kv.getValue()+"\t"+endOffset.get(kv.getKey())+"\n");
        }
        idxWriter.close();
    }
    private boolean containsTrueAltAllele(final List<VariantContext> VCs) {
        if ( VCs == null ) throw new IllegalArgumentException("The list of VariantContexts cannot be null");

        for ( final VariantContext vc : VCs ) {
            if ( vc.getNAlleles() > 2 )
                return true;
        }
        return false;
    }
    private VariantContext referenceBlockMerge(final List<VariantContext> VCs, final GenomeLocation loc,byte refAfter, final int end) {

        final VariantContext first = VCs.get(0);

        // ref allele and start
        final Allele refAllele;
        final int start;
        if ( loc == null || !loc.getContig().equals(first.getChr()) || first.getStart() >= loc.getStart() + 1) {
            start = first.getStart();
            refAllele = first.getReference();
        } else {
            start = loc.getStart() + 1;
            refAllele = Allele.create(refAfter, true);
        }
        // attributes
        final Map<String, Object> attrs = new HashMap<>(1);

        // genotypes
        final GenotypesContext genotypes = GenotypesContext.create();
        for ( final VariantContext vc : VCs ) {
            for ( final Genotype g : vc.getGenotypes() )
                genotypes.add(new GenotypeBuilder(g).alleles(GaeaGvcfVariantContextUtils.noCallAlleles(g.getPloidy())).make());
        }
        return new VariantContextBuilder("", first.getChr(), start, end, Arrays.asList(refAllele, VariantContextMerger.NON_REF_SYMBOLIC_ALLELE)).attributes(attrs).genotypes(genotypes).make();
    }
    private VCFHeader getVCFHeaderFromInput(Set<VCFHeader> headers) {
        Set<String> samples = getSampleList(headers);
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(headers, true);
        VCFHeader vcfHeader = new VCFHeader(headerLines, samples);

        headers.clear();
        samples.clear();
        headerLines.clear();

        return vcfHeader;
    }
    private Set<String> getSampleList(Set<VCFHeader> headers){
        Set<String> samples = new TreeSet<String>();
        for(VCFHeader header : headers){
            for ( String sample : header.getGenotypeSamples() ) {
                samples.add(GaeaGvcfVariantContextUtils.mergedSampleName(null, sample, false));
            }
        }
        return samples;
    }
    private void writeToDisk(String s, VariantContext vc) throws IOException {
        //写到哪里
        FileWriter w=new FileWriter((new File(vVcfPath)).getParent()+"/"+s+".mr2.lastVC");
        w.write(s);
        w.write("\n");
        w.write(vcfEncoder.encode(vc));
        w.write("\n");
        w.close();
        //最后在driver中要清理
    }
}
