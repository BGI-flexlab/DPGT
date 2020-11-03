package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.broadcast.Broadcast;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;

public class MergeToChrom implements Function2<Integer, Iterator<Integer>, Iterator<String>> {
//    private final String[] args;
    private final Integer cycle;
    private final DriverBC dBC;
    Logger log= LoggerFactory.getLogger(MergeToChrom.class);
    public MergeToChrom(Integer iter, Broadcast<DriverBC> dBC) {
        this.dBC = dBC.value();
        this.cycle=iter;
    }

    @Override public Iterator<String> call(Integer integer, Iterator<Integer> integerIterator) throws Exception {
        //write to outDir/chr*.vcf.gz

        ArrayList<String> r=new ArrayList<>();
        r.add("done");
        if(cycle.intValue()<=integer.intValue()){
            return r.iterator();
        }
        Logger logger = LoggerFactory.getLogger(MergeToChrom.class);

        String chr="";
        for(Map.Entry<String,Integer> entry:dBC.chrIndex.entrySet()){
            if(entry.getValue().equals(integer)){
                chr=entry.getKey();
                break;
            }
        }
        if(integer==25){
            chr="others";
        }
        File outFileName=new File(dBC.outputDir+"/"+chr+".vcf.gz");
        if(outFileName.exists()){
            return r.iterator();
        }
        BlockCompressedOutputStream out=new BlockCompressedOutputStream(dBC.outputDir+"/"+chr+".vcf.gz");
//        SeekableStream in = new SeekableFileStream(new File(dBC.outputDir+"/merged.vcfheader"));
//        VCFHeader header = VCFHeaderReader.readHeaderFrom(in);
        //write header
        String mergedHeaderFile=dBC.options.getOutDir()+"/merged.vcfheader";
        File hf=new File(mergedHeaderFile);
        if(hf.exists()){
            BufferedReader hfReader=new BufferedReader(new FileReader(mergedHeaderFile));
            String hfLine="";
            while((hfLine=hfReader.readLine())!=null){
                out.write((hfLine+"\n").getBytes());
            }
            hfReader.close();
        }else{
            logger.error("no header file");
            System.exit(1);
        }
        long start=1;
        long end=start+dBC.options.getRegionSize()-1;
        String firstFile="";
        ArrayList<String> justCopyList=new ArrayList<>();
        String lastFile="";
        for(int i=0;i<cycle;i++){
            File gtDir=new File(dBC.outputDir+"/genotype."+i);
            if(gtDir.exists() && gtDir.isDirectory()){
                boolean breakFlag=false;
                String[] gtLists=gtDir.list();
                String prefix="";
                for(int j=0;j<gtLists.length;j++){
                    if(gtLists[j].endsWith(".gz")){
                        String[] eles=gtLists[j].split("\\.");
                        if(eles.length==0){
                            logger.error("cannot found spot in file name:"+gtLists[0]);
                            System.exit(1);
                        }
                        prefix=eles[0];
                        break;
                    }
                }
                String[] poses=prefix.split("_");
                log.info("current pos:\t"+prefix);
                if(Long.valueOf(poses[1])<dBC.accumulateLength.get(integer)){
                    continue;
                }
                for(int j=0;j<dBC.options.getReducerNumber();j++){
                    File gtFile=new File(gtDir+"/"+prefix+"."+j+".gz");
                    if(gtFile.exists()){
                        BlockCompressedInputStream in=new BlockCompressedInputStream(gtFile);
                        String line1=in.readLine();
                        if(line1==null){
                            continue;
                        }
                        int tabNum=0;
                        String curChr="";
                        for(int k=0;k<line1.length();k++){
                            if(line1.charAt(k)=='\t'){
                                tabNum++;
                            }
                            if(tabNum==0){
                                curChr+=line1.charAt(k);
                            }else{
                                break;
                            }
                        }
                        if(dBC.chrIndex.get(curChr)==integer.intValue()){
                            justCopyList.add(gtFile.getAbsolutePath());
                        }else if(dBC.chrIndex.get(curChr)<integer.intValue()){
                            firstFile=gtFile.getAbsolutePath();
                        }else{
                            if(integer==25){
                                justCopyList.add(gtFile.getAbsolutePath());
                            }else {
                                lastFile = justCopyList.remove(justCopyList.size() - 1);
                                breakFlag = true;
                                in.close();
                                break;
                            }
                        }
                        in.close();
                    }
                }
                if(breakFlag){
                    break;
                }
            }
        }
        //process first file
        log.info("first file:\t"+firstFile);
        for(String file:justCopyList){
            log.info("middle file:\t"+file);
        }
        log.info("last file:\t"+lastFile);
        if(firstFile!=""){
            File gtFile=new File(firstFile);
            BlockCompressedInputStream in=new BlockCompressedInputStream(gtFile);
            String line1;
            while((line1=in.readLine())!=null) {
                int tabNum = 0;
                String curChr = "";
                for (int k = 0; k < line1.length(); k++) {
                    if (line1.charAt(k) == '\t') {
                        tabNum++;
                    }
                    if (tabNum == 0) {
                        curChr += line1.charAt(k);
                    } else {
                        break;
                    }
                }
                if (integer == 25) {
                    if(dBC.chrIndex.get(curChr)>=integer.intValue()){
                        out.write((line1 + "\n").getBytes());
                    }else{
                        continue;
                    }
                } else {
                    if (dBC.chrIndex.get(curChr) == integer.intValue()) {
                        out.write((line1 + "\n").getBytes());
                    } else if (dBC.chrIndex.get(curChr) < integer.intValue()) {
                        continue;
                    }else{
                        break;
                    }
                }
            }
            in.close();
        }
        //process copy files
        mergeBGZF(justCopyList,out);
        //process last file
        if(lastFile!=""){
            File gtFile=new File(lastFile);
            BlockCompressedInputStream in=new BlockCompressedInputStream(gtFile);
            String line1;
            while((line1=in.readLine())!=null) {
                int tabNum = 0;
                String curChr = "";
                for (int k = 0; k < line1.length(); k++) {
                    if (line1.charAt(k) == '\t') {
                        tabNum++;
                    }
                    if (tabNum == 0) {
                        curChr += line1.charAt(k);
                    } else {
                        break;
                    }
                }
                if (dBC.chrIndex.get(curChr) == integer.intValue()) {
                    out.write((line1 + "\n").getBytes());
                } else if (dBC.chrIndex.get(curChr) < integer.intValue()) {
                    continue;
                }else{
                    break;
                }
            }
            in.close();
        }
        out.close();
        //create tbi index
//        VCFCodec codec = new VCFCodec();
//        File chrVcf=new File(dBC.outputDir+"/"+chr+".vcf.gz");
//        TabixIndex tbiIndex= IndexFactory.createTabixIndex(new File(dBC.outputDir+"/"+chr+".vcf.gz"),codec, TabixFormat.VCF,null);
//        //Index tbiIndex=ctor.finalizeIndex(iterator.getPosition());
//
//        tbiIndex.writeBasedOnFeatureFile(chrVcf);
        return r.iterator();
    }
    void mergeBGZF(ArrayList<String> inputList, BlockCompressedOutputStream outFile) throws IOException {
        for(String file:inputList){
            BlockCompressedInputStream in=new BlockCompressedInputStream(new File(file));
            byte[] buf=new byte[10*1024*1024];
            int readSize=0;
            while((readSize=in.read(buf))>0){
                outFile.write(buf,0,readSize);
            }
            log.info("file merge done:\t"+file);
        }
        log.info("merge done");
    }
}
