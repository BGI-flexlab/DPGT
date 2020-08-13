package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

public class indexPartition extends org.apache.spark.Partitioner {
    @Override public int numPartitions() {
        return 0;
    }

    @Override public int getPartition(Object key) {
        return 0;
    }
}
