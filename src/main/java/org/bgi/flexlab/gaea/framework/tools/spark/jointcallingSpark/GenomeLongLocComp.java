package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;

import java.io.Serializable;
import java.util.Comparator;

public class GenomeLongLocComp
        implements Comparator<GenomeLongRegion>,Serializable {

    @Override public int compare(GenomeLongRegion o1, GenomeLongRegion o2) {

        if(o1.getStart()!=o2.getStart()){
            return (int) (o1.getStart()-o2.getStart());
        }else{
            return (int) (o1.getEnd()-o2.getEnd());
        }
    }
}
