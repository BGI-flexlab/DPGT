package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.broadcast.Broadcast;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class MergeToChrom implements Function2<Integer, Iterator<Integer>, Iterator<String>> {
//    private final String[] args;
    private final Integer cycle;
    private final DriverBC dBC;
    final Logger log= LoggerFactory.getLogger(MergeToChrom.class);
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
        String outputName=dBC.outputDir+"/"+chr+".vcf.gz";
        File outFileName=new File(outputName);
        File idxFileName=new File(outputName+".tbi");
        if(outFileName.exists() && idxFileName.exists()){
            return r.iterator();
        }
        if(!outFileName.exists()) {
            log.info("merging parallel files to",outFileName);
            BlockCompressedOutputStream out = new BlockCompressedOutputStream(outputName);
//        SeekableStream in = new SeekableFileStream(new File(dBC.outputDir+"/merged.vcfheader"));
//        VCFHeader header = VCFHeaderReader.readHeaderFrom(in);
            //write header
            String mergedHeaderFile = dBC.options.getOutDir() + "/merged.vcfheader";
            File hf = new File(mergedHeaderFile);
            if (hf.exists()) {
                BufferedReader hfReader = new BufferedReader(new FileReader(mergedHeaderFile));
                String hfLine;
                while ((hfLine = hfReader.readLine()) != null) {
                    out.write((hfLine + "\n").getBytes());
                }
                hfReader.close();
            } else {
                logger.error("no header file");
                System.exit(1);
            }
            long start = 1;
            long end = start + dBC.options.getRegionSize() - 1;
            String firstFile = "";
            ArrayList<String> justCopyList = new ArrayList<>();
            String lastFile = "";
            for (int i = 0; i < cycle; i++) {
                File gtDir = new File(dBC.outputDir + "/genotype." + i);
                if (gtDir.exists() && gtDir.isDirectory()) {
                    boolean breakFlag = false;
                    String[] gtLists = gtDir.list();
                    String prefix = "";
                    for (int j = 0; j < gtLists.length; j++) {
                        if (gtLists[j].endsWith(".gz")) {
                            String[] eles = gtLists[j].split("\\.");
                            if (eles.length == 0) {
                                logger.error("cannot found spot in file name:" + gtLists[0]);
                                System.exit(1);
                            }
                            prefix = eles[0];
                            break;
                        }
                    }
                    String[] poses = prefix.split("_");
                    if(poses.length<2){
                        continue;
                    }
                    log.info("current pos:\t" + prefix);
                    if (Long.valueOf(poses[1]) < dBC.accumulateLength.get(integer)) {
                        continue;
                    }
                    for (int j = 0; j < dBC.options.getReducerNumber(); j++) {
                        File gtFile = new File(gtDir + "/" + prefix + "." + j + ".gz");
                        if (gtFile.exists()) {
                            BlockCompressedInputStream in = new BlockCompressedInputStream(gtFile);
                            String line1 = in.readLine();
                            if (line1 == null) {
                                continue;
                            }
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
                                justCopyList.add(gtFile.getAbsolutePath());
                            } else if (dBC.chrIndex.get(curChr) < integer.intValue()) {
                                firstFile = gtFile.getAbsolutePath();
                            } else {
                                if (integer == 25) {
                                    justCopyList.add(gtFile.getAbsolutePath());
                                } else {
                                    if(justCopyList.size()>0) {
                                        lastFile = justCopyList.remove(justCopyList.size() - 1);
                                        breakFlag = true;
                                        in.close();
                                        break;
                                    }
                                    if(!firstFile.isEmpty()) {
                                        breakFlag = true;
                                        in.close();
                                        break;
                                    }
                                }
                            }
                            in.close();
                        }
                    }
                    if (breakFlag) {
                        break;
                    }
                }
            }
            //process first file
            log.info("first file:\t" + firstFile);
            for (String file : justCopyList) {
                log.info("middle file:\t" + file);
            }
            log.info("last file:\t" + lastFile);
            if (!firstFile.equals("")) {
                File gtFile = new File(firstFile);
                BlockCompressedInputStream in = new BlockCompressedInputStream(gtFile);
                String line1;
                while ((line1 = in.readLine()) != null) {
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
                        if (dBC.chrIndex.get(curChr) >= integer.intValue()) {
                            out.write((line1 + "\n").getBytes());
                        } else {
                            continue;
                        }
                    } else {
                        if (dBC.chrIndex.get(curChr) == integer.intValue()) {
                            out.write((line1 + "\n").getBytes());
                        } else if (dBC.chrIndex.get(curChr) < integer.intValue()) {
                            continue;
                        } else {
                            break;
                        }
                    }
                }
                in.close();
            }
            //process copy files
            mergeBGZF(justCopyList, out);
            //process last file
            if (!lastFile.equals("")) {
                File gtFile = new File(lastFile);
                BlockCompressedInputStream in = new BlockCompressedInputStream(gtFile);
                String line1;
                while ((line1 = in.readLine()) != null) {
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
                    } else {
                        break;
                    }
                }
                in.close();
            }
            out.close();
        }
        //create tbi index
//        VCFCodec codec = new VCFCodec();
//        File chrVcf=new File(dBC.outputDir+"/"+chr+".vcf.gz");
//        TabixIndex tbiIndex= IndexFactory.createTabixIndex(new File(dBC.outputDir+"/"+chr+".vcf.gz"),codec, TabixFormat.VCF,null);
//        //Index tbiIndex=ctor.finalizeIndex(iterator.getPosition());
//
//        tbiIndex.writeBasedOnFeatureFile(chrVcf);

        //create tbi index
        log.info("creating tbi index for:",outputName);
        VCFCodec codec = new VCFCodec();
        File output=new File(dBC.outputDir+"/"+chr+".vcf.gz");
        TabixIndexCreator ic=new TabixIndexCreator(TabixFormat.VCF);
//
        VariantContext lastContext = null;
        VariantContext currentContext;
        final Map<String, VariantContext> visitedChromos = new HashMap<String, VariantContext>();
        BlockCompressedInputStream bci=new BlockCompressedInputStream(output);
        AsciiLineReader lineReader = new AsciiLineReader(new BlockCompressedInputStream(output));
        AsciiLineReaderIterator iterator = new AsciiLineReaderIterator(lineReader);
        codec.readActualHeader(iterator);
        iterator.close();
        lineReader.close();
        String line;
        boolean readHeader=false;
        long position=0;
        while (true) {
            if(readHeader) {
                position = bci.getPosition();
            }
            if((line=bci.readLine())==null){
                break;
            }
            if(line.startsWith("#")){
                readHeader=true;
                continue;
            }
//            final long position2 = BlockCompressedFilePointerUtil.getBlockAddress(bci.getPosition());
            currentContext = codec.decode(line);
            checkSorted(lastContext, currentContext);
            //should only visit chromosomes once
            final String curChr = currentContext.getContig();
            final String lastChr = lastContext != null ? lastContext.getContig() : null;
            if(!curChr.equals(lastChr)){
                if(visitedChromos.containsKey(curChr)){
                    throw new RuntimeException("Input file must have contiguous chromosomes.");
                }else{
                    visitedChromos.put(curChr, currentContext);
                }
            }
            ic.addFeature(currentContext, position);
            lastContext = currentContext;
        }

//        iterator.close();
        Index idx=ic.finalizeIndex(bci.getPosition());
        idx.writeBasedOnFeatureFile(output);
        bci.close();
        return r.iterator();
    }


    private void checkSorted(VariantContext lastContext, VariantContext currentContext) {
        if (lastContext != null && currentContext.getStart() < lastContext.getStart() && lastContext.getChr().equals(currentContext.getChr()))
            throw new RuntimeException("Input file is not sorted by start position.");
    }

    void mergeBGZF(ArrayList<String> inputList, BlockCompressedOutputStream outFile) throws IOException {
        for(String file:inputList){
            BlockCompressedInputStream in=new BlockCompressedInputStream(new File(file));
            byte[] buf=new byte[10*1024*1024];
            int readSize;
            while((readSize=in.read(buf))>0){
                outFile.write(buf,0,readSize);
            }
            log.info("file merge done:\t"+file);
        }
        log.info("merge done");
    }
}
