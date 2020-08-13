package org.bgi.flexlab.gaea.tools.callsv;

import org.apache.hadoop.io.Writable;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapreduce.Partitioner;

/**
 * 重写ChrPartitionar类，指定每一个key所去的分区
 * @author Huifang Lu
 *
 */
public class ChrPartitionar extends Partitioner<NewMapKey, Writable>{

	public void configure(JobConf conf) { }

	/**
	 * 重写getPartition()方法，获取分区号，计算方法<br>
	 * 获取NewMapKey对象key中的chr，得到chr的hashCode值取绝对值，然后对numPartitions取余<br>
	 * Math.abs(key.getChr().hashCode()) % numPartitions;
	 * @param key NewMapKey类型的键
	 * @param value Writable类或者其子类对象
	 * @param numPartitions int类型，程序设置的reducer数目
	 * @return 分区号
	 * 
	 */
	public int getPartition(NewMapKey key, Writable value, int numPartitions) {		
		return Math.abs(key.getChr().hashCode()) % numPartitions;
	}

}
