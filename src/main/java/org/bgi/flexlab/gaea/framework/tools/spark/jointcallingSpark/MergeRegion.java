package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

public class MergeRegion implements FlatMapFunction<Iterator<GenomeLongRegion>, GenomeLongRegion>
         {
    public GenomeLongRegion gloc;
    public MergeRegion(GenomeLongRegion loc){
        gloc=loc;
    }
    @Override public Iterator<GenomeLongRegion> call(Iterator<GenomeLongRegion> genomeLocationIterator) throws Exception {

        Set<Long> realBreakpoints=new TreeSet<>();
        while(genomeLocationIterator.hasNext()){
            GenomeLongRegion loc=genomeLocationIterator.next();
            Long realEnd=loc.getEnd()>gloc.getEnd()?gloc.getEnd():loc.getEnd();
            for(Long i=loc.getStart();i<=realEnd;i++){
                realBreakpoints.add(i);
            }
        }
        long wStart=0;
        long wEnd=0;
        long lastPos=0;
        ArrayList<GenomeLongRegion> realRegion=new ArrayList<>();
        for(Long pos:realBreakpoints) {
            if(wStart==0) {
                wStart=pos;
                wEnd=pos;
            }else {
                if(pos-lastPos==1) {
                    wEnd=pos;
                }else {
                    GenomeLongRegion loc=new GenomeLongRegion(wStart,wEnd);
                    realRegion.add(loc);
                    wStart=pos;
                    wEnd=pos;
                }
            }
            lastPos=pos;
        }
        if(wStart!=0) {
            GenomeLongRegion loc=new GenomeLongRegion(wStart,wEnd);
            realRegion.add(loc);
        }
        realBreakpoints.clear();
        return realRegion.iterator();
    }
}
