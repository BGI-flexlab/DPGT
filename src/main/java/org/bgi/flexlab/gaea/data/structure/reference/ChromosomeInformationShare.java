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
package org.bgi.flexlab.gaea.data.structure.reference;

import org.bgi.flexlab.gaea.data.structure.memoryshare.BioMemoryShare;
import org.bgi.flexlab.gaea.util.SystemConfiguration;

/**
 * 染色体信息共享内存
 * 
 * @author ZhangYong
 *
 */
public class ChromosomeInformationShare extends BioMemoryShare {

	private int NonNbaselength = 0;

	public ChromosomeInformationShare() {
		super(Byte.SIZE / 4);
	}

	/**
	 * get base from reference position
	 * 
	 * @param pos
	 * @return base
	 */
	public byte getBinaryBase(int pos) {
		byte curr = getGA4GHBytes(pos, pos)[0];

		if ((pos & 0x1) == 0)
			return (byte) (curr & 0x0f);
		return (byte) ((curr >> 4) & 0x0f);
	}

	/**
	 * 获取碱基
	 */
	public char getBase(int pos) {
		return SystemConfiguration.getFastaAbb(getBinaryBase(pos));
	}

	public boolean isSNP(int pos) {
		byte posByte = getBinaryBase(pos);

		if (((posByte >> 3) & 0x1) == 0)
			return false;
		return true;
	}

	public boolean[] isSNPs(int start, int end) {
		if(end >= length)
			end = length - 1;
		
		byte[] bases = getGA4GHBytes(start, end);
		return isSNPs(bases, start, end);
	}
	
	public boolean[] isSNPs(byte[] bases,int start,int end){
		if(end >= length)
			end = length - 1;
		
		boolean[] snps = new boolean[end - start + 1];
		
		int deviation = (start & 0x1) == 1 ? 1 : 0;
		
		for(int i = start ; i < (end + 1) ; i++) {
			int posi = (i-start+deviation) / capacity;
			
			if((i & 0x1) == 0)
				snps[i-start] = ((bases[posi] >> 3) & 0x1) == 1;
			else
				snps[i-start] = ((bases[posi] >> 7) & 0x1) == 1;
		}
		
		return snps;
	}

	/**
	 * 获取染色体序列
	 * 
	 * @param start
	 *            从0开始
	 * @param end
	 * @return String 序列
	 */
	public String getGA4GHBaseSequence(int start, int end) {
		if(end >= length)
			end = length -1;
		
		byte[] bases = getGA4GHBytes(start, end);
		return getGA4GHBaseSequence(bases, start, end);
	}
	
	public String getGA4GHBaseSequence(byte[] bases,int start,int end){
		if(end >= length)
			end = length-1;
		
		StringBuilder seq = new StringBuilder();
		
		int deviation = (start & 0x1) == 1 ? 1 : 0;
			
		for(int i = start ; i < (end+1) ; i++){
			int posi = (i-start+deviation) / capacity;

			try {
				if ((i & 0x1) == 0)
					seq.append(SystemConfiguration.getFastaAbb(bases[posi] & 0x0f));
				else
					seq.append(SystemConfiguration.getFastaAbb((bases[posi] >> 4) & 0x0f));
			} catch(Exception e) {
				throw  new RuntimeException(chrName+"\tstart-end:" + start + "-" + end + "\ti:" + i + "\tposi:" +posi+"\t"+bases.length);
			}
		}
		
		return seq.toString();
	}
	
	public byte[] getGA4GHBaseBytes(int start) {
		return getGA4GHBaseSequence(start, start).getBytes();
	}

	public byte[] getGA4GHBaseBytes(int start, int end) {
		return getGA4GHBaseSequence(start, end).getBytes();
	}

	/**
	 * 转换byte数组为字符串
	 * 
	 * @param b
	 *            输入byte数组
	 * @return 返回相应的字符串
	 */
	public static String bytes2String(byte[] b) {
		StringBuilder sb = new StringBuilder(b.length);
		for (byte each : b) {
			sb.append((char) each);
		}
		return sb.toString();
	}

	/**
	 * byte转换为double类型
	 */
	public static double byteToDouble(byte[] b, int i) {
		long l = 0;
		for (int j = 0; j < 8; j++) {
			l |= (((long) (b[i + j] & 0xff)) << (8 * j));
		}
		return Double.longBitsToDouble(l);
	}

	public int getNonNbaselength() {
		if( NonNbaselength == 0) {
			for (int i = 0; i < length; i++) {
				if (getBase(i) != 'N') {
					NonNbaselength++;
				}
			}
		}
		return NonNbaselength;
	}
}
