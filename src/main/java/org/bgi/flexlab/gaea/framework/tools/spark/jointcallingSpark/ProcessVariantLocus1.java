package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.tribble.TabixFeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.spark.TaskContext;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.broadcast.Broadcast;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.BitSet;
import java.util.Iterator;


public class ProcessVariantLocus1 implements VoidFunction<Iterator<String>> {
    public GenomeLocation gloc;
    private final DriverBC dBC;
    public GenomeLongRegion region;
    public int iter;
    public int ｐarallelismNum;
    private static final Logger logger = LoggerFactory.getLogger(ProcessVariantLocus.class);

    public ProcessVariantLocus1(GenomeLocation loc, GenomeLongRegion region, int ｐarallelismNum, int iter, Broadcast<DriverBC> dBC) {
        this.gloc = loc;
        this.region = region;
        this.ｐarallelismNum = ｐarallelismNum;
        this.iter = iter;
        this.dBC = dBC.value();
    }

//    @Override
//    public void call(Iterator<String> stringIterator) throws Exception {
//        int id=TaskContext.getPartitionId();
//        Set<Long> realBreakpoints=new TreeSet<>();
//        //使用query方式试试
//        VCFCodec query_codec=new VCFCodec();
//
//        File rawInput=new File(dBC.options.getInputList());
//        String sortedInputGvcfList=dBC.options.getOutDir()+"/sorted."+rawInput.getName();
//        BufferedReader smallHeaderInfoReader=new BufferedReader(new FileReader(sortedInputGvcfList));
//        String headerInfoLine;
//        TabixFeatureReader sampleReader;
//
//        int start= gloc.getStart()+(gloc.getEnd() - gloc.getStart())/ｐarallelismNum*id;
//        int end=Math.min(gloc.getEnd(),gloc.getStart()+(gloc.getEnd() - gloc.getStart())/ｐarallelismNum*(id+1));
//        while((headerInfoLine=smallHeaderInfoReader.readLine())!=null){
//            sampleReader=new TabixFeatureReader(headerInfoLine,query_codec);
//            Iterator<VariantContext> it = sampleReader.query(gloc.getContig(), start, end);
//            while (it.hasNext()) {
//                final VariantContext vc = it.next();
//                if (vc.getNAlleles() > 2) {
//                    long curPosStart=dBC.accumulateLength.get(dBC.chrIndex.get(vc.getContig()))+vc.getStart();
//                    long curPosEnd=dBC.accumulateLength.get(dBC.chrIndex.get(vc.getContig()))+vc.getEnd();
////                    treeSet.add(new GenomeLongRegion(curPosStart,curPosEnd));
//                    long realEnd=region.getEnd()>curPosEnd?curPosEnd:region.getEnd();
//                    for(long i = curPosStart; i<=realEnd; i++){
//                        realBreakpoints.add(i);
//                    }
//                }
//            }
//        }
//        smallHeaderInfoReader.close();
//        String outputBP=dBC.options.getOutDir()+"/subbps."+ iter +"/part."+id;
//        File outOutFilePath=new File(outputBP);
//        if(outOutFilePath.exists()){
//            outOutFilePath.delete();
//        }
//        BufferedWriter bpsWriter=new BufferedWriter(new FileWriter(outputBP));
//
//        long wStart=0;
//        long wEnd=0;
//        long lastPos=0;
//        for(Long pos:realBreakpoints) {
//            if(wStart==0) {
//                wStart=pos;
//            }else {
//                if (pos - lastPos != 1) {
//                    bpsWriter.write(wStart+"\t"+wEnd+"\n");
//                    wStart = pos;
//                }
//            }
//            wEnd=pos;
//            lastPos=pos;
//        }
//        if(wStart!=0) {
//            bpsWriter.write(wStart+"\t"+wEnd+"\n");
//        }
//        bpsWriter.close();
//        realBreakpoints.clear();
//    }


    @Override
    public void call(Iterator<String> stringIterator) throws Exception {
        int id = TaskContext.getPartitionId();

        //使用query方式试试
        VCFCodec query_codec = new VCFCodec();

        File rawInput = new File(dBC.options.getInputList());
        String sortedInputGvcfList = dBC.options.getOutDir() + "/sorted." + rawInput.getName();
        BufferedReader smallHeaderInfoReader = new BufferedReader(new FileReader(sortedInputGvcfList));
        String headerInfoLine;
        TabixFeatureReader sampleReader;

        int start = gloc.getStart() + (gloc.getEnd() - gloc.getStart()) / ｐarallelismNum * id;
        int end = Math.min(gloc.getEnd(), gloc.getStart() + (gloc.getEnd() - gloc.getStart()) / ｐarallelismNum * (id + 1));

        BitSet bitSet = new BitSet(end - start + 1);
        long contigStart = dBC.accumulateLength.get(dBC.chrIndex.get(gloc.getContig()));

        while ((headerInfoLine = smallHeaderInfoReader.readLine()) != null) {
            sampleReader = new TabixFeatureReader(headerInfoLine, query_codec);
            Iterator<VariantContext> it = sampleReader.query(gloc.getContig(), start, end);
            while (it.hasNext()) {
                final VariantContext vc = it.next();
                if (vc.getNAlleles() > 2) {
                    int realEnd = region.getEnd() - contigStart > vc.getEnd() ? vc.getEnd() : (int) (region.getEnd() - contigStart);
                    for (int i = Math.max(0,vc.getStart() - start); i <= realEnd - start; i++) {
                        bitSet.set(i);
                    }
                }
            }
        }
        smallHeaderInfoReader.close();
        String outputBP = dBC.options.getOutDir() + "/subbps." + iter + "/part." + id;
        File outOutFilePath = new File(outputBP);
        if (outOutFilePath.exists()) {
            outOutFilePath.delete();
        }
        BufferedWriter bpsWriter = new BufferedWriter(new FileWriter(outputBP));

        int wStart = 0;
        int wEnd = 0;
        int lastPos = 0;
        long finalstart;
        long finalend;
        for (int pos = bitSet.nextSetBit(0); pos >= 0; pos = bitSet.nextSetBit(pos+1)) {
                if (wStart == 0) {
                    wStart = pos;
                } else {
                    if (pos - lastPos != 1) {
                         finalstart=contigStart + start+wStart;
                         finalend=contigStart + start+wEnd;
                        bpsWriter.write( finalstart+ "\t" + finalend + "\n");
                        wStart = pos;
                    }
                }
                wEnd = pos;
                lastPos = pos;
        }
        if (wStart != 0) {
             finalstart=contigStart + start+wStart;
             finalend=contigStart + start+wEnd;
            bpsWriter.write(finalstart + "\t" + finalend + "\n");
        }

        bpsWriter.close();
    }
}
