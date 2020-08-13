package org.bgi.flexlab.gaea.tools.callsv;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class ListComputer {
	
	/**
	 * 得到一个List的平均值 
	 * @param l 保存了insert size的List
	 * @return float型的平均值mean
	 */

	public static float getMean(List<Integer> l) {
		
		if(l.size() == 0) {
			return 0;
		}
		long sum = 0;
		for(Iterator<Integer> insertlist = l.iterator(); insertlist.hasNext();) {
			sum += insertlist.next();
		}
		
		return sum/l.size();
	}
	
	/**
	 * 获得一个列表的标准差
	 * @param l 保存了insert size的List
	 * @return float型的标准差
	 */
	
	public static float getStd(List<Integer> l) {
		float mean = getMean(l);
		return getUpLowStd(l, mean);
	}
	
	/**
	 * 获得Up或者Low列表的标准差
	 * @param l
	 * @param mean
	 * @return
	 */
	public static float getUpLowStd(List<Integer> l, float mean) {
		if(l.size() == 0) {
			return 0;
		}
		
		long sum_x2 = 0;
		for(Iterator<Integer> insertlist = l.iterator(); insertlist.hasNext();) {
			sum_x2 += Math.pow((insertlist.next() - mean), 2);
		}
		return (float) Math.sqrt(sum_x2/l.size());
	}
	
	/**
	 * 根据给定的平均值和标准差，删除列表中异常的（大于 mean+sd*sdTimes）的元素
	 * @param list 原始列表
	 * @param mean 平均值
	 * @param sd 标准差
	 * @param sdTimes 计算异常值时标准差的倍数
	 * @return 删除异常值后的列表
	 */
	public static List<Integer> delectOutlier(List<Integer> list, float mean, float sd, int sdTimes){
		List<Integer> l = new ArrayList<Integer>();
		for(Integer x : list) {
			if(x < mean+sd*sdTimes) {
				l.add(x);
			}
		}
		return l;
		
	}


}
