package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.tribble.TabixFeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import scala.Tuple2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;


public class ProcessVariantLocus implements PairFlatMapFunction<String, GenomeLongRegion, Integer> {
    private final ArrayList<GenomeLocation> regions=new ArrayList<>();
    private static final Logger logger = LoggerFactory.getLogger(ProcessVariantLocus.class);
    private final DriverBC dBC;
    public ProcessVariantLocus(ArrayList<GenomeLocation> regions, Broadcast<DriverBC> dBC) throws IOException {
        this.regions.addAll(regions);
        this.dBC=dBC.value();

    }

    @Override
    public Iterator<Tuple2<GenomeLongRegion,Integer>> call(String s) throws Exception {
        LinkedList<Tuple2<GenomeLongRegion,Integer>> variantsList = new LinkedList<>();
        if(s.startsWith("file://"))
            s=s.substring(7);

        //使用query方式试试
        VCFCodec query_codec=new VCFCodec();

        TabixFeatureReader sampleReader=new TabixFeatureReader(s,query_codec);
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
