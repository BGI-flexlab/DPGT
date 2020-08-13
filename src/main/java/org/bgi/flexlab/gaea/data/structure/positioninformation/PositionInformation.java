/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.structure.positioninformation;

import java.io.Serializable;
import java.util.Arrays;

public class PositionInformation implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -7761815558454832261L;
	
	/**
	 * 位点在参考基因组上的碱基编码
	 */
	private byte base;
	
	
	/**
	 * 唯一比对到该位点的所有read中的对应碱基（A/T/C/G）的各自数量之和数组
	 */
	private short[] uniReadBaseSum;
	
	/**
	 * 唯一比对到该位点的read中的对应碱基（A/T/C/G）的质量值之和数组
	 */
	private int[] qualityValueSum;
	
	/**
	 * 位点的测序深度
	 */
	private int seqDepth;
	
	/**
	 * 该位点至少被核苷酸唯一覆盖的次数
	 */
	private short uniCoveredCount;

	/**
	 * 该位点最佳比对到参考基因组上的总数
	 */
	private int repeatTime;

	/**
	 * 比对到该位点的所有read中的对应碱基（A/T/C/G）的各自数量之和数组
	 */
	private short[] allReadBaseSum;
	
	/**
	 * 经过过滤后的read的base sum
	 */
	private short[] readBaseSum;
	
	/**
	 * vcf 统计量
	 */
	private long baseQualSum;
	
	private int baseQualNum;
	
	private int[] sbstrandSum;
	
	private int[][] sameTypeQualityCount;

	/**
	 * 位点的坐标值
	 */
	private int position;

	/**
	 * 
	 */
	public PositionInformation() {
		allReadBaseSum = new short[4];
		uniReadBaseSum = new short[4];
		qualityValueSum = new int[4];
		base = (byte) 0xFF;
		
		baseQualSum = 0;
		baseQualNum = 0;
		sbstrandSum = new int[8];
		Arrays.fill(sbstrandSum, 0);
		
		readBaseSum = new short[4];
		sameTypeQualityCount = new int[4][95];
		for(int i = 0; i < 4; i++) {
			Arrays.fill(sameTypeQualityCount[i], 0);
		}
	}
	
	/*
	 * vcf 统计量
	 */
	public long getBaseQualSum() {
		return baseQualSum;
	}
	
	public void addBaseQualSum(int baseQual) {
		this.baseQualSum += baseQual; 
	}
	
	public int getBaseQualNum() {
		return baseQualNum;
	}
	
	public void addBaseQualNum(int num) {
		this.baseQualNum += num;
	}
	
	public int[] getSbstrandSum() {
		return sbstrandSum;
	}
	
	public void addSbstrandSum(int index, int sum) {
		sbstrandSum[index] += sum;
	}
	
	public int[][] getSameTypeQualityCount() {
		return sameTypeQualityCount;
	}
	
	public void addSameTypeQualityCount(int baseType, int qual, int sum) {
		sameTypeQualityCount[baseType][qual] += sum;
	}
	
	/**
	 * 获取该位点在参考基因组上的碱基编码
	 * @return
	 */
	public byte getBase() {
		return base;
	}

	/**
	 * 设置该位点位点在参考基因组上的碱基编码
	 * @param base
	 */
	public void setBase(byte base) {
		this.base = base;
	}
	
	/**
	 * 获取唯一比对到该位点的所有read中的对应碱基（A/T/C/G）的各自数量之和数组的指定元素
	 * @param index 数组下标
	 * @return
	 */
	public short getUniReadBaseSum(int index) {
		return uniReadBaseSum[index];
	}

	/**
	 * 设置唯一比对到该位点的所有read中的对应碱基（A/T/C/G）的各自数量之和数组的指定元素
	 * @param index
	 * @param value
	 */
	public void setUniReadBaseSum(int index, short value) {
		this.uniReadBaseSum[index] = value;
	}
	
	/**
	 * 获取唯一比对到该位点的read中的对应碱基（A/T/C/G）的质量值之和数组
	 * @return
	 */
	public int[] getQualityValueSum() {
		int[] q = new int[4];
		System.arraycopy(qualityValueSum, 0, q, 0, 4);
		return q;
	}

	/**
	 * 获取唯一比对到该位点的read中的对应碱基（A/T/C/G）的质量值之和数组的指定元素
	 * @param index 数组下标
	 * @return
	 */
	public int getQualityValueSum(int index) {
		return qualityValueSum[index];
	}

	/**
	 * 设置唯一比对到该位点的read中的对应碱基（A/T/C/G）的质量值之和数组的指定元素
	 * @param index 数组下标
	 * @param value
	 */
	public void setQualityValueSum(int index, int value) {
		this.qualityValueSum[index] = value;
	}

	/**
	 * 获取位点的测序深度
	 * @return
	 */
	public int getSeqDepth() {
		return seqDepth;
	}

	/**
	 * 设置位点的测序深度
	 * @param seqDepth
	 */
	public void setSeqDepth(int seqDepth) {
		this.seqDepth = seqDepth;
	}

	/**
	 * 获取该位点至少被核苷酸唯一覆盖的次数
	 * @return
	 */
	public short getUniCoveredCount() {
		return uniCoveredCount;
	}

	/**
	 * 设置该位点至少被核苷酸唯一覆盖的次数
	 * @param uniCoveredCount
	 */
	public void setUniCoveredCount(short uniCoveredCount) {
		this.uniCoveredCount = uniCoveredCount;
	}

	/**
	 * 获取该位点最佳比对到参考基因组上的总数
	 * @return
	 */
	public int getRepeatTime() {
		return repeatTime;
	}

	/**
	 * 设置该位点最佳比对到参考基因组上的总数
	 * @param repeatTime
	 */
	public void setRepeatTime(int repeatTime) {
		this.repeatTime = repeatTime;
	}

	/**
	 * 获取比对到该位点的所有read中的对应碱基（A/T/C/G）的各自数量之和数组的指定元素
	 * @param index 数组下标
	 * @return
	 */
	public short getAllReadBaseSum(int index) {
		return allReadBaseSum[index];
	}

	/**
	 * 获取比对到该位点的所有read中的对应碱基（A/T/C/G）的各自数量之和数组
	 * @return
	 */
	public short[] getAllReadBaseSum() {
		short[] a = new short[4];
		System.arraycopy(allReadBaseSum, 0, a, 0, 4);
		return a;
	}

	/**
	 * 设置比对到该位点的所有read中的对应碱基（A/T/C/G）的各自数量之和数组的指定元素
	 * @param index 数组下标
	 * @param value
	 */
	public void setAllReadBaseSum(int index, short value) {
		this.allReadBaseSum[index] = value;
	}

	/**
	 * 获取位点的坐标值
	 * @return
	 */
	public int getPosition() {
		return position;
	}

	/**
	 * 设置位点的坐标值       
	 * @param position
	 */
	public void setPosition(int position) {
		this.position = position;
	}

	/**
	 * 
	 * @param index
	 * @return
	 */
	public short getReadBaseSum(int index) {
		return readBaseSum[index];
	}

	/**
	 * 
	 * @param index
	 * @param value
	 */
	public void setReadBaseSum(int index, short value) {
		this.readBaseSum[index] = value;
	}
	
}
