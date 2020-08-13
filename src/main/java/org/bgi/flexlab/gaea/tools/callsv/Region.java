package org.bgi.flexlab.gaea.tools.callsv;

import java.util.ArrayList;
import java.util.List;

/**
 * Region类，将划分好的区域封装，这个类中包含了描述region的成员<br>
 * 包括： 区域编号regId，所属的染色体chr，区域的起始和终止位置regStart和regEnd<br>
 *       这个区域的reads数目，碱基数目，以及比对到这个区域的正向和反向reads数目，和readsId列表
 * @author Huifang Lu
 *
 */
public class Region {
	
	/**
	 * 区域编号
	 */
	private int regId;
	/**
	 * 区域所属的染色体
	 */
	private String chr;
	/**
	 * 区域的起始位置
	 */
	private int regStart;
	/**
	 * 区域的终止位置
	 */
	private int regEnd;
	/**
	 * 区域中上一条reads的起始位置；
	 */
	private int pStart;
	/**
	 * 比对到区域的reads数目
	 */
	private int regReadNum;
	/**
	 * 所有reads的碱基总数
	 */
	private int baseNum;
	/**
	 * 正向比对的reads数目
	 */
	private int positiveOriNum;
	/**
	 * 反向比对的reads数目
	 */
	private int negativeOriNum;
	/**
	 * 所有readsID的列表
	 */
	private List<String> regReads;
	
	/**
	 * 初始化一个Region对象，成员变量都是默认值
	 */
	public Region(){
		this.regId = 0;
		this.regEnd = 0;
		this.chr = null;
		this.regStart = 1000000000;
		this.pStart = 0;
		this.regReadNum = 0;
		this.baseNum = 0;
		this.positiveOriNum = 0;
		this.negativeOriNum = 0;
		this.regReads = new ArrayList<String>();
	}
	
	/**
	 * 根据参数初始化Region对象，成员变量用参数赋值
	 * @param regId 区域编号
	 * @param chr 区域所属的染色体
	 * @param regStart 区域起始位置
	 * @param regEnd 区域终止位置
	 * @param regReadNum 区域的reads数
	 * @param positiveOriNum 正向比对的reads数
	 * @param negativeOriNum 反向比对的reads数
	 * @param regReads 区域readsID的列表
	 */
	public Region(int regId, String chr, int regStart, int regEnd, int pStart, int regReadNum, int positiveOriNum, int negativeOriNum, List<String> regReads){
		this.regId = regId;
		this.chr = chr;
		this.regStart = regStart;
		this.regEnd = regEnd;
		this.pStart = pStart;
		this.regReadNum = regReadNum;
		this.positiveOriNum = positiveOriNum;
		this.negativeOriNum = negativeOriNum;
		this.regReads = regReads;
	}
	
	/**
	 * 使用参数初始化Region对象，chr，regEnd, regStart由参数f指定，其它为默认值
	 * @param f 一条reads比对信息
	 */
	public Region(SamWritable f){
		this.regId = 0;
		this.chr = f.getChr();
		this.regEnd = f.getEnd();
		this.regStart = f.getStart();
		this.pStart = f.getStart();
		this.regReadNum = 0;
		this.baseNum = 0;
		this.positiveOriNum = 0;
		this.negativeOriNum = 0;
		this.regReads = new ArrayList<String>();
	}
	
	public int getRegId() {
		return regId;
	}
	
	public void setRegId(int regId) {
		this.regId = regId;
	}
	
	public String getChr() {
		return chr;
	}
	
	public void setChr(String chr) {
		this.chr = chr;
	}
	
	public int getRegStart() {
		return regStart;
	}
	
	public void setRegStart(int regStart) {
		this.regStart = regStart;
	}
	
	public int getRegEnd() {
		return regEnd;
	}
	
	public void setRegEnd(int regEnd) {
		this.regEnd = regEnd;
	}
	

	public int getpStart() {
		return pStart;
	}

	public void setpStart(int pStart) {
		this.pStart = pStart;
	}

	public int getRegReadNum() {
		return regReadNum;
	}
	
	public void setRegReadNum(int regReadNum) {
		this.regReadNum = regReadNum;
	}
	
	public int getBaseNum() {
		return baseNum;
	}

	public void setBaseNum(int baseNum) {
		this.baseNum = baseNum;
	}

	public int getPositiveOriNum() {
		return positiveOriNum;
	}
	
	public void setPositiveOriNum(int positiveOriNum) {
		this.positiveOriNum = positiveOriNum;
	}
	
	public int getNegativeOriNum() {
		return negativeOriNum;
	}
	
	public void setNegativeOriNum(int negativeOriNum) {
		this.negativeOriNum = negativeOriNum;
	}
	
	public List<String> getRegReads() {
		return regReads;
	}
	
	public void setRegReads(List<String> regReads) {
		this.regReads = regReads;
	}
	
	public int getRegLength() {
		int length = regEnd - regStart + 1;
		return length;
	}
	
	public float getRegCoverage() {
		int length = getRegLength();
		if(length == 0)
			return 0;
		else
			return baseNum/length;
	}

	@Override
	public int hashCode() {
		return chr.hashCode() + regId*31;
	}

	@Override
	public boolean equals(Object obj) {
		if(!(obj instanceof Region))
			throw new ClassCastException("obj can not cast to Region class!");
		Region r = (Region)obj;
		return this.chr.equals(r.chr) && this.regId==r.regId;
	}

	@Override
	public String toString() {
		return regId + "\t" + chr + "\t" + regStart + "\t" + regEnd + "\t" + pStart + "\t" + regReadNum + "\t" + positiveOriNum + "+" + negativeOriNum + "-";
	}

	/**
	 * 根据reads的比对信息更新Region对象的成员变量值
	 * @param r 一条reads的比对信息
	 */
	public void updateReg(SamWritable r) {
		int start = Math.min(this.regStart, r.getStart());
		this.setRegStart(start);
		
		int end = Math.max(this.regEnd, r.getEnd());
		this.setRegEnd(end);
		
		this.regReadNum ++;
		this.baseNum = this.baseNum + r.getReadLen();
		this.pStart = r.getStart();
		
		if(r.getStrand().equals("+"))
			this.positiveOriNum ++;
		else
			this.negativeOriNum ++;
		
		
		regReads.add(r.getReadName());
		
	}
	
	
	public String firstToString() {
		return chr + "\t" + regEnd + "\t" + positiveOriNum + "+" + negativeOriNum + "-";
	}
	
	public String secondToString() {
		return chr + "\t" + regStart + "\t" + positiveOriNum + "+" + negativeOriNum + "-";
	}

}
