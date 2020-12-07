package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.vcf.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.util.LongAccumulator;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalWriter;
import org.bgi.flexlab.gaea.util.Utils;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;

@SuppressWarnings({"ALL",
        "ResultOfMethodCallIgnored"})
public class JointCallingSpark {
    private final static String INPUT_ORDER = "input.name.order";
    private final static String INPUT_LIST = "input.gvcf.list";// added by gc
    private final static String Window_File = "window.file.path";// window file path
    private static final Logger logger = LoggerFactory.getLogger(JointCallingSpark.class);
    private static final LinkedHashMap<String, Integer> chrIndex=new LinkedHashMap<>();
//    public static ArrayList<BufferedReader> sampleReaders=new ArrayList<>();
//    public static ArrayList<BufferedReader> sampleReadersForMR2=new ArrayList<>();//HashMap可以深拷贝，不担心互相受影响
    private static GenomeLocation parseRegionFromString(String targetRegion) {
//        GenomeLocation gloc;
        String ele="";
        ArrayList<String> eles=new ArrayList<>();
        for(int i=0;i<targetRegion.length();i++){
            if(Character.isLetterOrDigit(targetRegion.charAt(i))){
                ele+=targetRegion.charAt(i);
            }else{
                eles.add(ele);
                ele="";
                if(eles.size()==3){
                    break;
                }
            }
        }
        if(eles.size()==2){
            if(ele!=""){
                eles.add(ele);
            }
        }
        if(eles.size()!=3){
            logger.error("fail to parse targetRegion, please check the options");
        }
        String chr=eles.get(0);
        if(chrIndex.containsKey(chr)){
            return new GenomeLocation(chr,chrIndex.get(chr),Integer.parseInt(eles.get(1)),Integer.parseInt(eles.get(2)));
        }else if(chrIndex.containsValue(Integer.parseInt(chr))){
            for(Map.Entry<String,Integer> kv:chrIndex.entrySet()){
                if(kv.getValue()==Integer.parseInt(chr)){
                    return new GenomeLocation(kv.getKey(),kv.getValue(),Integer.parseInt(eles.get(1)),Integer.parseInt(eles.get(2)));
                }
            }
            logger.error("code error");
        }else{
            logger.error("no such chromosome name/id:"+chr+", please check the options");
        }
        logger.error("fail to parse targetRegion, please check the options");
        return null;
    }
    public static void main(String[] args) throws IOException {
        String HEADER_DEFAULT_PATH = "vcfheader";
        String MERGER_HEADER_INFO = "vcfheaderinfo";
        ArrayList<ArrayList<String>> multiMapSampleNames=new ArrayList<>();

        JointCallingSparkOptions options = new JointCallingSparkOptions();
        options.parse(args);
        String vVcfPath;
        LinkedHashMap<String,Integer> sampleIndex=new LinkedHashMap<>();
        Map<String, String> pathSample = new HashMap<>();


        VCFHeader header = null;
        LinkedHashMap<Integer, String> contigs = null;
        SparkConf conf = new SparkConf().setAppName("JointCallingSpark");

        conf.set("spark.serializer", "org.apache.spark.serializer.JavaSerializer");
        conf.set("spark.rdd.compress","true");
        JavaSparkContext sc = new JavaSparkContext(conf);   //打开spark环境
        Configuration hadoopConf=sc.hadoopConfiguration();
        File tmpDir=new File(options.getOutDir());
        if(!tmpDir.exists()){
            tmpDir.mkdirs();
        }


        Configuration hadoop_conf=new Configuration(hadoopConf);
//        MultipleVCFHeaderForJointCalling multiVcfHeader = new MultipleVCFHeaderForJointCalling();
        logger.warn("before get Header");
//        ArrayList<Path> pathList=new ArrayList<>();
        VCFHeader mergeHeader=null;
//        for(String a:options.getInputStringList()){
//            pathList.add(new Path(a));
//        }
        File headerDir=new File(options.getOutDir()+"/headers");


        File vcfheaderFile=new File(options.getOutDir()+"/vcfheader");
        File vcfheaderInfoFile=new File(options.getOutDir()+"/vcfheaderinfo");
        if(!vcfheaderFile.exists() || !vcfheaderInfoFile.exists() || !headerDir.exists()) {
            if(!headerDir.exists()){
                headerDir.mkdirs();
            }else{
                if(!headerDir.isDirectory()){
                    headerDir.delete();
                    headerDir.mkdirs();
                }
            }
            String inputList=options.getInputList();
            if(!inputList.startsWith("file://")){
                inputList="file://"+inputList;
            }
            JavaRDD<String> gvcfSamples = sc.textFile(inputList, options.getMapperNumber());
            gvcfSamples.mapPartitionsWithIndex(new ProcessHeader(options), true).collect();
        }
        //create vcfheader and vcfheaderinfo by small files in directory "headers"
        BufferedWriter vcfInfoWriter=new BufferedWriter(new FileWriter(options.getOutDir()+"/"+MERGER_HEADER_INFO));
        int inputIndex=0;
        LinkedHashSet<VCFHeader> headers=new LinkedHashSet<>();
        TreeSet<String> sampleNames=new TreeSet<>();
        hadoop_conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getOutDir() + "/"+HEADER_DEFAULT_PATH);
        hadoop_conf.set(MERGER_HEADER_INFO, options.getOutDir()+"/"+MERGER_HEADER_INFO);
        HashMap<String,String> sampleNamePath=new HashMap<>();
        HashMap<String,String> sampleNameRPath=new HashMap<>();
        for(int i=0;i<options.getMapperNumber();i++){
            File smallHeader=new File(options.getOutDir()+"/headers/vcfheader"+i);
            File smallHeaderInfo=new File(options.getOutDir()+"/headers/vcfPathName"+i);
            if(!smallHeader.exists()){
                logger.error("file not exists,"+smallHeader.getAbsolutePath());
                System.exit(1);
            }
            if(!smallHeaderInfo.exists()){
                logger.error("file not exists,"+smallHeaderInfo.getAbsolutePath());
                System.exit(1);
            }

            BufferedReader smallHeaderInfoReader=new BufferedReader(new FileReader(smallHeaderInfo.getAbsolutePath()));
            String headerInfoLine;
            while((headerInfoLine=smallHeaderInfoReader.readLine())!=null){
                String[] eles=headerInfoLine.split("\t");
                sampleNamePath.put(eles[2],eles[0]);
                sampleNameRPath.put(eles[2],eles[1]);
            }
            smallHeaderInfoReader.close();
            VCFLocalLoader vcfLL=new VCFLocalLoader(smallHeader.getAbsolutePath());
            headers.add(vcfLL.getHeader());
            sampleNames.addAll(vcfLL.getHeader().getSampleNamesInOrder());
            if(inputIndex>=1000){
                VCFHeader vcfHeader=new VCFHeader(VCFUtils.smartMergeHeaders(headers, true),sampleNames);
                headers.clear();
                headers.add(vcfHeader);
            }
        }
        if (headers.size() >= 1) {
            mergeHeader=new VCFHeader(VCFUtils.smartMergeHeaders(headers, true),sampleNames);
            headers.clear();
        }
        inputIndex=0;
        File rawInput=new File(options.getInputList());
        String sortedInputGvcfList=options.getOutDir()+"/sorted."+rawInput.getName();
        BufferedWriter sortedInputWriter=new BufferedWriter(new FileWriter(sortedInputGvcfList));
        for(String sampleName:mergeHeader.getSampleNamesInOrder()){
            if(!sampleNamePath.containsKey(sampleName)){
                logger.error("code error");
                System.exit(1);
            }
            sortedInputWriter.write(sampleNamePath.get(sampleName)+"\n");
            vcfInfoWriter.write(sampleNameRPath.get(sampleName)+"\t"+inputIndex+"\t"+sampleName+"\n");
            sampleIndex.put(sampleName,inputIndex);
            pathSample.put(sampleNameRPath.get(sampleName), sampleName);
            inputIndex++;
        }
        sortedInputWriter.close();
        vcfInfoWriter.close();
        VCFLocalWriter mergedVCFHeaderWriter=new VCFLocalWriter(options.getOutDir()+"/"+HEADER_DEFAULT_PATH,false,false);
        mergedVCFHeaderWriter.writeHeader(mergeHeader);
        mergedVCFHeaderWriter.close();
        if(!sortedInputGvcfList.startsWith("file://")) {
            sortedInputGvcfList = "file://" + sortedInputGvcfList;
        }
//        options.setInputList(sortedInputGvcfList);
        JavaRDD<String> sortedGvcfSamples = sc.textFile(sortedInputGvcfList, options.getMapperNumber());
        FileWriter virtualVCF=new FileWriter(new File(options.getOutDir()+"/virtual.vcf"));
        String fileList=sortedInputGvcfList;
        if(fileList.startsWith("file://")){
            fileList=fileList.substring(7);
        }
        BufferedReader headerReader=new BufferedReader(new FileReader(options.getOutDir()+"/vcfheader"));

        String tmpline;
        while ((tmpline = headerReader.readLine()) != null) {
            if (tmpline.startsWith("#CHROM")) {
                String writeLine = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVirtualSample\n";
                virtualVCF.write(writeLine);
                break;
            } else {
                virtualVCF.write(tmpline);
                virtualVCF.write("\n");
            }
        }
        virtualVCF.close();
        //使用虚拟vcf文件，提取contig信息
        SeekableStream in=new SeekableFileStream(new File(options.getOutDir()+"/virtual.vcf"));
        header = VCFHeaderReader.readHeaderFrom(in);
        in.close();
        if (header == null)
            throw new RuntimeException("header is null !!!");
        contigs = new LinkedHashMap<>();
        for (VCFContigHeaderLine line : header.getContigLines()) {
            contigs.put(line.getContigIndex(), line.getID());
            chrIndex.put(line.getID(),line.getContigIndex());
        }

        //创建窗口文件
        String win_out_file = options.getOutDir() + "/windows.bed";
        // Path raw_win_file=new Path(options.getWinFile());//get windows file path from
        // command
        if (win_out_file.startsWith("file://")) {
            win_out_file = win_out_file.substring(7);
        }
        File winOutFile = new File(win_out_file);
        if (!winOutFile.exists()) {
            winOutFile.createNewFile();
        }
        FileWriter win_out = new FileWriter(winOutFile); // output stream ready for write
        int window_size = options.getWindowsSize();
        for (Map.Entry<Integer, String> entry : contigs.entrySet()) {
            String chr = entry.getValue();
            int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
            int start = 1;
            int end = -1;
            while (true) {
                end = start + window_size - 1;
                if (end > contigLength) {
                    end = contigLength;
                }
                String write_line = chr + "\t" + start + "\t" + end + "\n";
                win_out.write(write_line);
                start += window_size;
                if (start > contigLength) {
                    break;
                }
            }
        }
        win_out.close();
        HashMap<String,String> confMap=new HashMap<>();
        confMap.put(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getOutDir() + "/"+HEADER_DEFAULT_PATH);
        confMap.put(MERGER_HEADER_INFO, options.getOutDir()+"/"+MERGER_HEADER_INFO);
        ArrayList<String> inputListArr = new ArrayList<>();
        BufferedReader gvcfList=new BufferedReader(new FileReader(fileList));
        String tmpLine;
        while ((tmpLine = gvcfList.readLine()) != null) {
            String[] eles = tmpLine.split("/");
            String rName = eles[eles.length - 1];
            String sampleName = pathSample.get(rName);
            inputListArr.add(sampleName);
        }
        gvcfList.close();
        String[] inputListStringArr = new String[inputListArr.size()];
        for (int i = 0; i < inputListArr.size(); i++) {
            inputListStringArr[i] = inputListArr.get(i);
        }
        hadoop_conf.set(INPUT_ORDER, Utils.join(",", inputListStringArr));
        confMap.put(INPUT_ORDER, Utils.join(",", inputListStringArr));

        inputListArr.clear();
        hadoop_conf.set(INPUT_LIST, sortedInputGvcfList);
        confMap.put(INPUT_LIST, sortedInputGvcfList);
        // create windows bed file based on window size
        hadoop_conf.set(Window_File, options.getOutDir() + "/windows.bed");
        confMap.put(Window_File, options.getOutDir() + "/windows.bed");
        vVcfPath=options.getOutDir()+"/virtual.vcf";
        logger.warn("after get Header");

        int[] index=new int[options.getMapperNumber()];
        for(int i=0;i<options.getMapperNumber();i++){
            index[i]=i;
        }
        for(List<String> samples:sortedGvcfSamples.collectPartitions(index)){
            ArrayList<String> indexSamples=new ArrayList<>();
            for(String sample:samples){
                String[] eles = sample.split("/");
                String rName = eles[eles.length - 1];
                String sampleName = pathSample.get(rName);
                indexSamples.add(sampleName);
            }
            multiMapSampleNames.add(indexSamples);
        }
        int iter=0;

//        Broadcast<ArrayList<BufferedReader>> sampleReadersForMR2BC=sc.broadcast(sampleReadersForMR2);
        HashMap<Integer,Long> accumulateLength=new HashMap<>();
        long totalLength=0;
        for(int ii=0;ii<chrIndex.size();ii++){
            long chrLength=header.getSequenceDictionary().getSequence(ii).getSequenceLength();
            totalLength+=chrLength;
            accumulateLength.put(ii,totalLength-chrLength);
        }
        accumulateLength.put(chrIndex.size(),totalLength);
        SeekableStream in2=new SeekableFileStream(new File(options.getOutDir()+"/virtual.vcf"));
        VCFHeader virtualHeader = VCFHeaderReader.readHeaderFrom(in2);
        Set<VCFHeaderLine> gvcfHeaderMetaInfo=virtualHeader.getMetaDataInInputOrder();
        VCFHeaderVersion version=null;
        for (final VCFHeaderLine line : gvcfHeaderMetaInfo) {
            if (VCFHeaderVersion.isFormatString(line.getKey())) {
                version = VCFHeaderVersion.toHeaderVersion(line.getValue());
                break;
            }
        }

        Broadcast<DriverBC> dBC=sc.broadcast(new DriverBC(options.getOutDir(),multiMapSampleNames,chrIndex,options, vVcfPath, sampleIndex, pathSample,accumulateLength,virtualHeader,version));
        Long step=options.getRegionSize()-1;
        Long longStart=1L;
        Long longEnd=longStart+step;
        GenomeLongRegion region=new GenomeLongRegion(longStart,longEnd);
        Long cycleEnd=totalLength;
        if(options.getTargetRegion()!=null){
            GenomeLocation targetRegion=parseRegionFromString(options.getTargetRegion());
            region.setStart(accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getStart());
            if(accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getEnd()-accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getStart()>step){
                region.setEnd(region.getStart()+step);
            }else{
                region.setEnd(accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getEnd());
            }
            if(accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getEnd()<cycleEnd){
                cycleEnd=accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getEnd();
            }
        }
        final LongAccumulator totalVariantsNum=sc.sc().longAccumulator("Variants num");
        while(true){
            //结束条件1：基因组区域处理结束
            if(region.getStart()>cycleEnd){
                break;
            }
            ArrayList<GenomeLocation> regions=new ArrayList<>();
            long startLoc=0;
            long endLoc=0;
            int startChr=-1;
            int endChr=-1;
            for(int ii=0;ii<accumulateLength.size();ii++){
                if(accumulateLength.get(ii)>region.getStart()){
                    if(startChr==-1) {
                        startChr=ii-1;
                        startLoc = region.getStart() - accumulateLength.get(ii - 1);
                    }
                }
                if(accumulateLength.get(ii)>=region.getEnd()){
                    if(endChr==-1) {
                        endChr=ii-1;
                        endLoc = region.getEnd() - accumulateLength.get(ii - 1);
                        break;
                    }
                }
            }
            if(region.getEnd()>=accumulateLength.get(accumulateLength.size()-1)){
                endChr=chrIndex.size()-1;
                endLoc=header.getSequenceDictionary().getSequence(endChr).getSequenceLength();
            }
            if(endChr==startChr){
                regions.add(new GenomeLocation(contigs.get(startChr),(int)startLoc,(int)endLoc));
            }else{
                for(int ii=startChr;ii<=endChr;ii++){
                    if(ii==startChr){
                        regions.add(new GenomeLocation(contigs.get(ii),(int)startLoc,header.getSequenceDictionary().getSequence(ii).getSequenceLength()));
                    }else if(ii==endChr){
                        regions.add(new GenomeLocation(contigs.get(ii),1,(int)endLoc));
                    }else{
                        regions.add(new GenomeLocation(contigs.get(ii),1,header.getSequenceDictionary().getSequence(ii).getSequenceLength()));
                    }
                }
            }
            //提取变异位置
            logger.info("region:\t"+region);
            logger.info("regions size:\t"+regions.size());
            for(GenomeLocation curRegion:regions){
                logger.info(curRegion.toString());
                if(curRegion.getContigIndex()>25){
                    break;
                }
            }
            JavaPairRDD<GenomeLongRegion,Integer> variantsRegion = sortedGvcfSamples.flatMapToPair(new ProcessVariantLocus(region,regions,dBC));
            //分区
            String outputBP=options.getOutDir()+"/bps."+String.valueOf(iter);
            File bpDir=new File(outputBP);
            if(!bpDir.exists() || !bpDir.isDirectory()) {

                JavaPairRDD<GenomeLongRegion, Integer> partitionedRegion = variantsRegion.partitionBy(new GenomeLocPartitioner(options.getMapperNumber(), region));
                JavaRDD<GenomeLongRegion> mergedRegion = partitionedRegion.keys().mapPartitions(new MergeRegion(region)).sortBy(x -> x.getStart(), true, 1);
                mergedRegion.saveAsTextFile("file://" + outputBP);
            }
            int codeInspectorIter=0;
            if(bpDir.isDirectory()){
                String[] files=bpDir.list();
                for(String file:files){
                    if(file.startsWith("part")){
                        codeInspectorIter++;
                    }
                }
            }
            if(codeInspectorIter>1){
                logger.error("code error, only one BPs file should be generated");
                System.exit(1);
            }
            String mergedAllBPs="";
            if(bpDir.isDirectory()){
                String[] files=bpDir.list();
                for(String file:files){
                    if(file.startsWith("part")){
                        mergedAllBPs=bpDir.getAbsolutePath()+"/"+file;
                    }
                }
            }
            BufferedReader bp_reader=new BufferedReader(new FileReader(mergedAllBPs));
            int bpIter=0;
            String bpRegion=bp_reader.readLine();
            bpIter++;
            if(bpRegion==null){
                bpIter++;
                region.setStart(region.getEnd()+1);
                region.setEnd(region.getStart()+step);
                continue;
            }else{
                while(bp_reader.readLine() !=null){
                    bpIter++;
                }
            }
            bp_reader.close();
            int variantsNumInEachReducer=bpIter/options.getReducerNumber();
            bp_reader=new BufferedReader(new FileReader(mergedAllBPs));
            //bpPartition: record end position of each genotype container process region
            ArrayList<Long> bpPartition=new ArrayList<>();
            bpIter=0;
            String lastBpRegion=null;
            while((bpRegion=bp_reader.readLine())!=null){
                bpIter++;
                if(bpIter%variantsNumInEachReducer==0 && bpPartition.size()<options.getReducerNumber()-1){
                    String[] eles=bpRegion.split("\t");
                    bpPartition.add(Long.parseLong(eles[0]));
                }
                lastBpRegion=bpRegion;
            }
            bp_reader.close();
            String[] bpEles=lastBpRegion.split("\t");
            bpPartition.add(Long.parseLong(bpEles[0])+1);
            if(bpPartition.size()!=options.getReducerNumber()){
                System.out.println("code error: bp partition error");
                System.exit(1);
            }
            String combineOutLocal=options.getOutDir()+"/combine."+String.valueOf(iter);
            String genotypeOut=options.getOutDir()+"/genotype."+String.valueOf(iter);
            File genotypeOutPath=new File(genotypeOut);
            boolean doGenotype=true;
            if(!genotypeOutPath.exists()){
                genotypeOutPath.mkdirs();
            }else{
                if(!genotypeOutPath.isDirectory()){
                    genotypeOutPath.delete();
                }else{
                    doGenotype=false;
                }
            }

            File combineOutLocalFile=new File(combineOutLocal);
            if(!combineOutLocalFile.exists()){
                combineOutLocalFile.mkdirs();
                if(doGenotype) {
                    sortedGvcfSamples.foreachPartition(new CombineVariants(region, regions, mergedAllBPs, confMap, dBC, iter, bpPartition));
                }
            }else{
                if(!combineOutLocalFile.isDirectory()){
                    combineOutLocalFile.delete();
                    combineOutLocalFile.mkdirs();
                    if(doGenotype) {
                        sortedGvcfSamples.foreachPartition(new CombineVariants(region, regions, mergedAllBPs, confMap, dBC, iter, bpPartition));
                    }
                }
            }
            ArrayList<String> tmpStr=new ArrayList<>();
            for(int i=0;i<options.getReducerNumber()*3;i++){
                tmpStr.add("abc"+String.valueOf(i));
            }
            JavaRDD<String> nonsenseRDD=sc.parallelize(tmpStr,options.getReducerNumber());
            if(doGenotype) {
                JavaRDD<String> variantsNum = nonsenseRDD.mapPartitionsWithIndex(new ParseCombineAndGenotypeGVCFs(region, regions, args, mergedAllBPs, confMap, dBC, iter, bpPartition,totalVariantsNum), true);
                variantsNum.collect();
            }

            iter++;
            region.setStart(region.getEnd()+1);
            region.setEnd(region.getStart()+step);
            if(!options.isKeepCombine()) {
                if (combineOutLocalFile.isDirectory()) {
                    String[] eles = combineOutLocalFile.list();
                    for (String ele : eles) {
                        File eleFile = new File(combineOutLocalFile.getAbsolutePath() + "/" + ele);
                        eleFile.delete();
                    }
                    combineOutLocalFile.delete();
                }
            }
            logger.info("current total variants:\t"+totalVariantsNum.value());
        }
        //create files of each contains whole chromosome data
        if(options.isMergeChrom()){
            //merge chr1-22,X,Y,M as default
            ArrayList<Integer> tmpList=new ArrayList<>();
            for(int i=1;i<=100;i++){
                tmpList.add(i);
            }
            JavaRDD<Integer> nonsenseRDD2=sc.parallelize(tmpList,26);
            nonsenseRDD2.mapPartitionsWithIndex(new MergeToChrom(iter,dBC),true).collect();
        }
        sc.stop();
    }

}
