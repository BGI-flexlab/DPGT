package org.bgi.flexlab.gaea.tools.callsv;

import org.apache.hadoop.io.WritableComparable;
import org.apache.hadoop.io.WritableComparator;

/**
 * Reducer的key区分标准类，只用chr来做区分
 * @author Huifang Lu
 *
 */
public class ReduceGroupingComparator extends WritableComparator{
	
	/**
	 * 构造方法
	 */
	public ReduceGroupingComparator(){
		super(NewMapKey.class, true);
	}

	@SuppressWarnings("rawtypes")
	@Override
	public int compare(WritableComparable a, WritableComparable b) {
		NewMapKey o1 = (NewMapKey) a;
		NewMapKey o2 = (NewMapKey) b;
	
		return o1.getChr().compareTo(o2.getChr());
	}
	

}
