package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang.RandomStringUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.bgi.flexlab.gaea.data.structure.dbsnp.DbsnpShare;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.data.structure.variant.VariantContextMerger;
import org.bgi.flexlab.gaea.tools.jointcalling.VariantAnnotatorEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.ReferenceConfidenceVariantContextMerger;
import org.bgi.flexlab.gaea.util.Utils;
import org.seqdoop.hadoop_bam.LazyVCFGenotypesContext;
import org.seqdoop.hadoop_bam.VariantContextWithHeader;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import scala.Int;
import scala.Tuple2;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.zip.GZIPInputStream;


public class ExtractVariants implements PairFlatMapFunction<Iterator<String>, Integer, Integer> {
    public GenomeLocation gloc;
    public String bpPath;
    public String vVcfPath;
    //    public Set<VCFHeaderLine> gvcfHeaderMetaInfo;
//    public VCFHeaderVersion version = null;
    public JointCallingSparkOptions options = new JointCallingSparkOptions();
    public final static String INPUT_LIST = "input.gvcf.list";
    HashMap<String, String> confMap=new HashMap<>();
    public String outputDir;
    Logger logger;
    public VCFEncoder vcfEncoder=null;
    public DriverBC dBC;
    public int iter;
    //    public static VCFHeader mergedHeader=null;
    public ExtractVariants(GenomeLocation region, String outputBP, HashMap<String, String> confMap,
                           Broadcast<DriverBC> dBC,int iter) throws
            IOException {
        gloc=region;
        bpPath=outputBP;
        this.dBC=dBC.value();
        this.options=dBC.value().options;
        this.confMap=confMap;
        this.outputDir=dBC.value().outputDir;
        this.confMap.put(INPUT_LIST, options.getInputList());
        logger = LoggerFactory.getLogger(ExtractVariants.class);
        this.vVcfPath=dBC.value().vVcfPath;
        this.iter=iter;

    }
    public class VcComp implements Comparator<VariantContext>,Serializable{

        @Override public int compare(VariantContext o1, VariantContext o2) {
            return o1.getStart()-o2.getStart();
        }
    }
    Comparator<VariantContext> comparator4 = new VcComp();
    @Override public Iterator<Tuple2<Integer, Integer>> call(Iterator<String> StringIterator) throws Exception {
        //与map的setup功能类似
        Configuration conf=new Configuration();
        Path path = new Path(outputDir+"/vcfheader");
        SeekableStream in2 = WrapSeekable.openPath(path.getFileSystem(conf), path);
        VCFHeader mergedHeader = VCFHeaderReader.readHeaderFrom(in2);
        if(mergedHeader==null){
            logger.error("header is null, "+in2.getSource());
            System.exit(1);
        }
        if(vcfEncoder==null)
            vcfEncoder=new VCFEncoder(mergedHeader, true, true);
        Configuration hadoop_conf=new Configuration();
        for(String k:confMap.keySet()){
            hadoop_conf.set(k,confMap.get(k));
        }
//        ArrayList<Tuple2<GenomeLocation,VariantContext>> vcList=new ArrayList<>();
        ArrayList<Tuple2<Integer,Integer>> vcList=new ArrayList<>();
        Set<String> mapSamplePath=new HashSet<>();
        Set<VCFHeaderLine> gvcfHeaderMetaInfo=null;
        HashMap<String, Integer> sampleIndex=new HashMap();
        String sampleIndexFile=options.getOutDir()+"/vcfheaderinfo";
        Path sampleIndexFilePath=new Path(sampleIndexFile);
        BufferedReader indexReader = new BufferedReader(new InputStreamReader(sampleIndexFilePath.getFileSystem(hadoop_conf).open(sampleIndexFilePath)));
        String indexLine;
        Logger logger= LoggerFactory.getLogger(ExtractVariants.class);
        Map<String,String> pathSample=new HashMap();
        ArrayList<String> sampleNames=new ArrayList<String>();
        Integer mapSMtagInt=0;
        int totalSampleSize=0;
        VCFHeaderVersion version = null;
        HashMap<Integer,Set<String>> mapIndexSample=new HashMap<>();
        ArrayList<VCFCodec> codec = new ArrayList<VCFCodec>();
        ArrayList<String> mapGvcfList=new ArrayList<String>();
        while(StringIterator.hasNext()) {
            String gvcfPath = StringIterator.next();
            mapGvcfList.add(gvcfPath);
        }
        for(String gvcfPath:mapGvcfList) {
            Path path2=new Path(gvcfPath);
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
        indexReader.close();
        HashMap<String,Integer> chrIndexs = null;
        VCFHeader merged_header = null;
        VCFHeader mapMergedHeader=null;
        String mapSMtag=null;
        HashMap<Integer, String> contigs = new HashMap<Integer, String>();
        Map<String, Integer> contigDict = new HashMap<String, Integer>();
        GenomeLocationParser parser = null;
        ReferenceShare genomeShare = null;
        DbsnpShare dbsnpShare = null;
        LazyVCFGenotypesContext.HeaderDataCache vcfHeaderDataCache =
                new LazyVCFGenotypesContext.HeaderDataCache();
        SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
        VariantAnnotatorEngine annotationEngine;
        List<String> annotationsToUse = new ArrayList<>();

        List<String> annotationGroupsToUse = new ArrayList<>(
                Arrays.asList(new String[] { StandardAnnotation.class.getSimpleName() }));

        BufferedReader[] gvcfBufReader=new BufferedReader[mapGvcfList.size()];
        ArrayList<VariantContext> curSamplesVC=new ArrayList<VariantContext>();

        Integer sISize=sampleIndex.size();
        int ii=0;
        ArrayList<TabixFeatureReader> samplesReader=new ArrayList<>();
        SeekableStream in=new SeekableFileStream(new File(options.getTmpOut()+"/virtual.vcf"));
        merged_header = VCFHeaderReader.readHeaderFrom(in);//从header文件中读取header
        in.close();
        for(String gvcfPath:mapGvcfList){
            String[] eles=gvcfPath.split("/");
            mapSamplePath.add(eles[eles.length-1]);
            Path path2=new Path(gvcfPath);
            if(!pathSample.containsKey(eles[eles.length-1])) {
                logger.error("no such path in vcfHeaderInfo");
            }
            String sampleName=pathSample.get(eles[eles.length-1]);
            sampleNames.add(sampleName);
            mapSMtagInt+=sampleIndex.get(sampleName);
            gvcfHeaderMetaInfo=merged_header.getMetaDataInInputOrder();
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
            mapIndexSample.put(ii,curSample);
            tmp_codec.setVCFHeader(new VCFHeader(gvcfHeaderMetaInfo,curSample), version);
            codec.add(tmp_codec);
            logger.warn(ii+"\tcurrent index");
            if(gvcfPath.startsWith("file://")) {
                gvcfPath=gvcfPath.substring(7);
            }
            BufferedReader reader=null;
//            VCFCodec query_codec=new VCFCodec();
            TabixFeatureReader sampleReader=new TabixFeatureReader(gvcfPath,tmp_codec);
            Iterator<VariantContext> it=sampleReader.query(gloc.getContig(),gloc.getStart(),gloc.getEnd());
            samplesReader.add(sampleReader);


        }
        gvcfHeaderMetaInfo=merged_header.getMetaDataInInputOrder();

        mapSamplePath.clear();


        logger.warn("sampleNames Size:\t"+sampleNames.size());
        logger.warn("mapGvcfList Size:\t"+mapGvcfList.size());
        pathSample.clear();
        sampleIndex.clear();
        mapSMtagInt+=totalSampleSize;
        mapMergedHeader=new VCFHeader(gvcfHeaderMetaInfo,sampleNames);
        vcfHeaderDataCache.setHeader(mapMergedHeader);

        mapSMtag= Utils.join("_",mapMergedHeader.getSampleNamesInOrder());
        if(mapMergedHeader.getSampleNamesInOrder().size()==1) {
            mapSMtag=mapSMtag+"_";
        }
        if(merged_header == null)
            throw new RuntimeException("header is null !!!");
        logger.warn("after get merged header");
        contigs = new HashMap<>();
        chrIndexs = new HashMap<>();
        for (VCFContigHeaderLine line : merged_header.getContigLines()) {
            chrIndexs.put(line.getID(), line.getContigIndex());
        }
        for (VCFContigHeaderLine line : merged_header.getContigLines()) {
            contigs.put(line.getContigIndex(), line.getID());
        }

        parser = new GenomeLocationParser(merged_header.getSequenceDictionary());

        logger.warn("after parser");

        genomeShare = new ReferenceShare();
        genomeShare.loadChromosomeList(options.getReference());
        dbsnpShare = new DbsnpShare(options.getDBSnp(), options.getReference());
        dbsnpShare.loadChromosomeList(options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
        for (final VCFContigHeaderLine contig : merged_header.getContigLines())
            contigDict.put(contig.getID(), contig.getContigIndex());
        logger.warn("after contigDict Map done");
        annotationEngine = new VariantAnnotatorEngine(annotationGroupsToUse, annotationsToUse,
                Collections.<String>emptyList());
        annotationEngine.initializeDBs(options.getDBSnp() != null);
        logger.warn("setup done");

        //与map的map功能类似
        String winFilePath=hadoop_conf.get(dBC.Window_File);
        if(winFilePath.startsWith("file://")) {
            winFilePath=winFilePath.substring(7);
        }
        String win_line;
        Path firstPath=new Path(mapGvcfList.get(0));
        FileSystem gvcf_fs=firstPath.getFileSystem(hadoop_conf);
        for(int i=0;i!=mapGvcfList.size();i++) {
            String gvcfFile=mapGvcfList.get(i);
            if(gvcfFile.startsWith("file://")) {
                gvcfFile=gvcfFile.substring(7);
            }
            BufferedReader reader=null;
            if(gvcfFile.endsWith(".gz")) {
                if(gvcfFile.startsWith("/user")) {
                    reader=new BufferedReader(new InputStreamReader(new GZIPInputStream(gvcf_fs.open(new Path(gvcfFile)))));
                }else {
                    reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gvcfFile))));
                }
            }else {
                reader = new BufferedReader(new InputStreamReader(new FileInputStream(gvcfFile)));
            }
            gvcfBufReader[i]=reader;
        }
        //method 3, merge code from MapperHuge
        BufferedReader win_reader=new BufferedReader(new FileReader(winFilePath));
        String last_contig=null;
        ChromosomeInformationShare ref=null;
        Path pBPpath=new Path(bpPath);
        BufferedReader bp_reader=new BufferedReader(new InputStreamReader(pBPpath.getFileSystem(hadoop_conf).open(pBPpath)));
        String winLine=bp_reader.readLine();
        String[] bpeles=winLine.split(":");
        Integer bpcontig=chrIndexs.get(bpeles[0]);
        Integer bpstart,bpend;
        if(bpeles[1].contains("-")){
            String[] eles2=bpeles[1].split("-");
            bpstart=Integer.valueOf(eles2[0]);
            bpend=Integer.valueOf(eles2[1]);
        }else{
            bpstart=bpend=Integer.valueOf(bpeles[1]);
        }
        ArrayList<Integer> bpStarts=new ArrayList<Integer>();
        ArrayList<Integer> bpEnds=new ArrayList<Integer>();
        Boolean bpFileEnd=false;
        Boolean logOut=false;
        Integer lastLogPoint=0;
        Integer last_end=0;
        Integer last_pos=0;
        Integer last_bp_start=-1;
        Integer last_bp_end=-1;
        Set<Integer> last_end_breakpoints=new TreeSet<Integer>();
        ArrayList<Iterator<VariantContext>> sampleVCsIter=new ArrayList<>();

        String contig=gloc.getContig();
        Integer start=gloc.getStart();
        Integer end=gloc.getEnd()>start+dBC.options.getWindowsSize()-1?start+dBC.options.getWindowsSize()-1:gloc.getEnd();
        ArrayList<Iterator<VariantContext>> samplesIt=new ArrayList<>();
        ArrayList<ArrayList<VariantContext>> samplesVcList=new ArrayList<>();
        int[] VcIndex=new int[mapGvcfList.size()];
        HashMap<Integer,VariantContext> lastVC=new HashMap<>();
        HashMap<String,Long> samplesOffset=new HashMap<>();
        while(end<=gloc.getEnd()) {
            Integer curLogPoint=(int)start/100000;
            if(!curLogPoint.equals(lastLogPoint)) {
                logOut=true;
            }else {
                logOut=false;
            }
            lastLogPoint=curLogPoint;
            if(logOut)
                System.out.println(formatter.format(new Date())+"\tprocess region\t"+contig+"\t"+start+"\t"+end);
            String chr = contig;

            Integer winChrID=chrIndexs.get(contig);
            if(last_bp_start!=-1){
                bpStarts.add(start);
                bpEnds.add(last_bp_end);
                last_bp_start=-1;
                last_bp_end=-1;
            }
            //将目标处理区域内的bp信息提取出来
            if(logOut)
                System.out.println(formatter.format(new Date())+"\tbp start\t"+contig+"\t"+start+"\t"+end);
            while(!bpFileEnd) {
                if(bpcontig.equals(winChrID)) {
                    if(bpstart<=end && bpend>=start) {
                        bpStarts.add(bpstart);
                        bpEnds.add(bpend);
                        winLine=bp_reader.readLine();
                        if(winLine==null) {
                            bpFileEnd=true;
                            bp_reader.close();
                            break;
                        }
                        bpeles=winLine.split(":");
                        bpcontig=chrIndexs.get(bpeles[0]);
                        if(bpeles[1].contains("-")){
                            String[] eles2=bpeles[1].split("-");
                            bpstart=Integer.valueOf(eles2[0]);
                            bpend=Integer.valueOf(eles2[1]);
                        }else{
                            bpstart=bpend=Integer.valueOf(bpeles[1]);
                        }
                    }else if(bpstart>end) {
                        if(bpend>end){
                            last_bp_start=bpstart;
                            last_bp_end=bpend;
                        }
                        break;
                    }else {
                        winLine=bp_reader.readLine();
                        if(winLine==null) {
                            bpFileEnd=true;
                            bp_reader.close();
                            break;
                        }
                        bpeles=winLine.split(":");
                        bpcontig=chrIndexs.get(bpeles[0]);
                        if(bpeles[1].contains("-")){
                            String[] eles2=bpeles[1].split("-");
                            bpstart=Integer.valueOf(eles2[0]);
                            bpend=Integer.valueOf(eles2[1]);
                        }else{
                            bpstart=bpend=Integer.valueOf(bpeles[1]);
                        }
                    }
                }else if(bpcontig>winChrID) {
                    break;
                }else {
                    winLine=bp_reader.readLine();
                    if(winLine==null) {
                        bpFileEnd=true;
                        bp_reader.close();
                        break;
                    }
                    bpeles=winLine.split(":");
                    bpcontig=chrIndexs.get(bpeles[0]);
                    if(bpeles[1].contains("-")){
                        String[] eles2=bpeles[1].split("-");
                        bpstart=Integer.valueOf(eles2[0]);
                        bpend=Integer.valueOf(eles2[1]);
                    }else{
                        bpstart=bpend=Integer.valueOf(bpeles[1]);
                    }
                }
            }
            if(logOut)
                System.out.println(formatter.format(new Date())+"\tbp end\t"+contig+"\t"+start+"\t"+end);
            LinkedList<VariantContext> region_vcs=new LinkedList<VariantContext>();
            Set<Integer> end_breakpoints=new TreeSet<Integer>();
            if(logOut) {
                System.out.println(formatter.format(new Date())+"\tbefore query");
            }
            for(int i=0;i!=curSamplesVC.size();i++) {
                VariantContext tmpVC=curSamplesVC.get(i);
                if(tmpVC.getChr().equals(chr) && tmpVC.getStart()<=end && tmpVC.getEnd()>=start) {
                    region_vcs.add(tmpVC);
                    if(tmpVC.getNAlleles()>2) {
                        for(int tmpPos=tmpVC.getStart();tmpPos<=tmpVC.getEnd();tmpPos++) {
                            end_breakpoints.add(tmpPos);
                        }
                    }else {
                        end_breakpoints.add(tmpVC.getEnd());
                    }
                    if(tmpVC.getEnd()<=end) {
                        curSamplesVC.remove(i);
                        i--;
                    }
                }
            }
            for(int i=0;i!=mapGvcfList.size();i++) {
                String line=null;
                String samplePath=mapGvcfList.get(i);

                String sampleName = sampleNames.get(i);
                Integer sampleIdx=dBC.sampleIndex.get(sampleName);
                long offset=0;
                while((line=gvcfBufReader[i].readLine())!=null) {
                    if(line.startsWith("#")) {
                        offset+=line.length()+1;
                        continue;
                    }
                    final VariantContext v = codec.get(i).decode(line);
                    CommonInfo info = v.getCommonInfo();
                    if(!info.hasAttribute("SM"))
                        info.putAttribute("SM",sampleName);
                    if(v.getChr().equals(chr)) {
                        if(v.getStart()<=end && v.getEnd()>=start) {
                            if(v.getNAlleles()>2) {
                                for(int tmpPos=v.getStart();tmpPos<=v.getEnd();tmpPos++) {
                                    end_breakpoints.add(tmpPos);
                                }
                            }else {
                                end_breakpoints.add(v.getEnd());
                            }
                            region_vcs.add(v);
                            if(v.getEnd()>end) {
                                curSamplesVC.add(v);
                                break;
                            }else{
                                offset+=line.length()+1;
                            }
                        }else if(v.getStart()>end) {
                            curSamplesVC.add(v);
                            break;
                        }else if(v.getEnd()<start) {
                            offset+=line.length()+1;
                            continue;
                        }else {
                        }
                    }else if(contigDict.get(v.getChr()) > contigDict.get(chr) ){
                        curSamplesVC.add(v);
                        break;
                    }else{
                        offset+=line.length()+1;
                        continue;
                    }
                }
                if(iter>0) {
                    samplesOffset.put(sampleName, samplesOffset.get(sampleName) + offset);
                }else{
                    samplesOffset.put(sampleName,offset);
                }
            }
            if(logOut) {
                System.out.println(formatter.format(new Date())+"\tafter query");
            }
            Collections.sort(region_vcs,comparator4);
            //do combineGVCFS

            if(!contig.equals(last_contig)) {
                ref=genomeShare.getChromosomeInfo(contig);
                last_pos=0;
            }else {
                last_pos=start;
            }
            ArrayList<Integer> realBreakpointStart=new ArrayList<Integer>();
            ArrayList<Integer> realBreakpointEnd=new ArrayList<Integer>();
            int bpIndex=0;
            for(Integer pos:end_breakpoints) {
                Boolean whether_keep=false;
                for(;bpIndex!=bpStarts.size();bpIndex++) {
                    if(last_pos<=bpEnds.get(bpIndex) && pos>=bpStarts.get(bpIndex)) {
                        whether_keep=true;
                        break;
                    }else if(last_pos>bpEnds.get(bpIndex)) {
                        continue;
                    }else if(pos<bpStarts.get(bpIndex)) {
                        break;
                    }
                }
                if(whether_keep) {
                    realBreakpointStart.add(last_pos);
                    realBreakpointEnd.add(pos);
                    if(pos>end){
                        last_end_breakpoints.add(pos);
                    }
                }
                last_pos=pos+1;
            }

            if(region_vcs.size()>0) {
                if(logOut) {
                    System.out.println(formatter.format(new Date())+"\tvc number in this region:\t"+region_vcs.size()+"\tcurSamplesVC:\t"+curSamplesVC.size()+"\tbpStarts:\t"+bpStarts.size()+"\tend_breakpoints:\t"+end_breakpoints.size()+"\trealBreakpoints:\t"+realBreakpointStart.size());
                }
            }
            end_breakpoints.clear();
            bpStarts.clear();
            bpEnds.clear();
            if(!contig.equals(last_contig)) {
                last_end=0;
            }
            //last_pos=0;
            int tmpIter=0;
            if(logOut) {
                System.out.println(formatter.format(new Date()) + "\tbefore merge");
            }
            for(int i=0;i!=realBreakpointStart.size();i++) {
                Integer posStart=realBreakpointStart.get(i);
                Integer pos=realBreakpointEnd.get(i);

                if(pos>end) {
                    break;
                }
                if(pos<start) {
                    continue;
                }
                List<VariantContext> stoppedVCs = new ArrayList<>();
                for(int j=0;j<region_vcs.size();j++) {
                    final VariantContext vc = region_vcs.get(j);
                    //the VC for the previous state will be stopped if its position is previous to the current position or it we've moved to a new contig
                    if ( vc.getStart() <= pos) {

                        if(vc.getEnd()<posStart) {
                            region_vcs.remove(j);
                            j--;
                        }else {
                            stoppedVCs.add(vc);
                        }
                        // if it was ending anyways, then remove it from the future state
                        if ( vc.getEnd() == pos) {
                            //region_vcws.samples.removeAll(vc.getSampleNames());
                            if(j>=0) {
                                region_vcs.remove(j);
                                j--;
                            }
                        }
                    }else {
                        break;
                    }
                }
                if ( !stoppedVCs.isEmpty()) {

                    final GenomeLocation gLoc =parser.createGenomeLocation(chr, last_end);
                    final VariantContext mergedVC;
                    final Byte refBase =(byte) ref.getBase(gLoc.getStart());
                    int curStart=last_end+1;
                    if ( containsTrueAltAllele(stoppedVCs) ) {
                        mergedVC = ReferenceConfidenceVariantContextMerger.mapMerge(stoppedVCs, parser.createGenomeLocation(chr, posStart), (byte) ref.getBase(posStart-1), false, false, annotationEngine);
                        tmpIter++;

                    }else {
                        mergedVC = referenceBlockMerge(stoppedVCs, gLoc,refBase, pos);
                    }
                    if(mergedVC==null) {
                        continue;
                    }

                    CommonInfo info = mergedVC.getCommonInfo();
                    if(!info.hasAttribute("SM"))
                        info.putAttribute("SM",mapSMtagInt);
                    HashMap<String, Object> maps = new HashMap<>();
                    maps.putAll(info.getAttributes());
                    info.setAttributes(maps);
                    GenomeLocation vcLoc=new GenomeLocation(mergedVC.getContig(),chrIndexs.get(mergedVC.getContig()),mergedVC.getStart(),mergedVC.getEnd());
                    vcList.add(new Tuple2<>(mergedVC.getStart(),1));
                    stoppedVCs.clear();
                }
                last_end=pos;
            }
            if(logOut) {
                System.out.println(formatter.format(new Date()) + "\tafter merge");
            }
            last_contig=contig;
            end_breakpoints.clear();
            start=end+1;
            end=gloc.getEnd()>end+dBC.options.getWindowsSize()?end+dBC.options.getWindowsSize():gloc.getEnd();
            if(start>end){
                break;
            }
        }
        String filename= RandomStringUtils.randomAlphanumeric(6);
        File offsetFile=new File(filename);
        while(offsetFile.exists()){
            filename= RandomStringUtils.randomAlphanumeric(6);
            offsetFile=new File(filename);
        }
        FileWriter offsetWrite=new FileWriter(offsetFile);
        for(Map.Entry<String,Long> entry:samplesOffset.entrySet()){
            offsetWrite.write(entry.getKey()+"\t"+entry.getValue());
        }
        offsetWrite.close();
        win_reader.close();

        return vcList.iterator();
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
    private VCFHeader getVCFHeaderFromInput(Set<VCFHeader> headers) throws IOException{
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
