package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;

import java.io.Serializable;
import java.util.Comparator;

public class GenomeLocComp
        implements Comparator<GenomeLocation>, Serializable {

    @Override public int compare(GenomeLocation o1, GenomeLocation o2) {
        if(o1.getContigIndex()!=o2.getContigIndex()){
            System.out.println("code error");
            System.exit(88);
        }
        if(o1.getStart()!=o2.getStart()){
            return o1.getStart()-o2.getStart();
        }else{
            return o1.getEnd()-o2.getEnd();
        }
    }
}
