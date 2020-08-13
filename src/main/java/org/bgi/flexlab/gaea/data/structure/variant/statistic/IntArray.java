package org.bgi.flexlab.gaea.data.structure.variant.statistic;

import java.util.Arrays;

public class IntArray {
	private int[] array;

	public IntArray(int length,int value){
		array = new int[length];
		Arrays.fill(array, 0);
	}
	
	public void incr(VariantEnum data,int add){
		array[data.ordinal()] += add;
	}
	
	public void incr(VariantEnum data){
		array[data.ordinal()]++;
	}
	
	public int get(VariantEnum data){
		return get(data.ordinal());
	}
	
	public int get(int index){
		return array[index];
	}
	
	public int length(){
		return array.length;
	}
	
	public int[] get(){
		return array;
	}
}
