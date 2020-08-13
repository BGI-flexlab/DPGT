package org.bgi.flexlab.gaea.tools.mapreduce.jointcallingprepare;

import org.apache.hadoop.mapreduce.Partitioner;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;

import java.util.Map;

public class WindowsBasedTestPartitioner<T> extends Partitioner<WindowsBasedWritable, T> {

	@Override
	public int getPartition(WindowsBasedWritable key, T v, int numPartitioner) {
		//chr1 1M-2M
		Integer pos=key.getPosition().get();
		Long curPos=pos+JointCallingPrepareMapper.accumulate.get(key.getChromosomeIndex());
		Long band=(long)(JointCallingPrepareMapper.refLength*1.1/numPartitioner);
		if(band==0){
			System.out.println("input data is empty or un-recognized");
			System.exit(-1);
		}
		return (int)(curPos / band);
	}
}