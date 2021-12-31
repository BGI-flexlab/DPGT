package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import org.apache.spark.Partitioner;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;

public class GenomeLocPartitioner extends Partitioner {
    public int partitions;
    public GenomeLongRegion gloc;
    public long partitionRange;
    public GenomeLocPartitioner(int parts,GenomeLongRegion region){
        this.partitions=parts;
        this.gloc=region;
        partitionRange=(gloc.getEnd()-gloc.getStart())/partitions;
    }
    @Override
    public int numPartitions() {
        return partitions;
    }
    /*
    todo
     */
    @Override
    public int getPartition(Object key) {
        GenomeLongRegion glocKey=(GenomeLongRegion)key;

        if(glocKey.getStart()-gloc.getStart()<0){
            return 0;
        }else if(glocKey.getStart()>gloc.getEnd()){
            return partitions-1;
        }else {
            int part=(int)((glocKey.getStart() - gloc.getStart()) / partitionRange);
            if(part>=partitions){
                return partitions-1;
            }else{
                return part;
            }
        }
    }
}
