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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.util;

import java.util.HashMap;
import java.util.Map;

public class SystemConfiguration {
	/**
	 * BAM flags
	 */
	/*
	 * ! @abstract the read is paired in sequencing, no matter whether it is
	 * mapped in a pair
	 */
	public static final int BAM_FPAIRED = 1;
	/* ! @abstract the read is mapped in a proper pair */
	public static final int BAM_FPROPER_PAIR = 2;
	/*
	 * ! @abstract the read itself is unmapped; conflictive with
	 * BAM_FPROPER_PAIR
	 */
	public static final int BAM_FUNMAP = 4;
	/* ! @abstract the mate is unmapped */
	public static final int BAM_FMUNMAP = 8;
	/* ! @abstract the read is mapped to the reverse strand */
	public static final int BAM_FREVERSE = 16;
	/* ! @abstract the mate is mapped to the reverse strand */
	public static final int BAM_FMREVERSE = 32;
	/* ! @abstract this is read1 */
	public static final int BAM_FREAD1 = 64;
	/* ! @abstract this is read2 */
	public static final int BAM_FREAD2 = 128;
	/* ! @abstract not primary alignment */
	public static final int BAM_FSECONDARY = 256;
	/* ! @abstract QC failure */
	public static final int BAM_FQCFAIL = 512;
	/* ! @abstract optical or PCR duplicate */
	public static final int BAM_FDUP = 1024;

	/*
	 * CIGAR operations.
	 */
	/* ! @abstract CIGAR: M = match or mismatch */
	public static final int BAM_CMATCH = 0;
	/* ! @abstract CIGAR: I = insertion to the reference */
	public static final int BAM_CINS = 1;
	/* ! @abstract CIGAR: D = deletion from the reference */
	public static final int BAM_CDEL = 2;
	/* ! @abstract CIGAR: N = skip on the reference (e.g. spliced alignment) */
	public static final int BAM_CREF_SKIP = 3;
	/*
	 * ! @abstract CIGAR: S = clip on the read with clipped sequence present in
	 * qseq
	 */
	public static final int BAM_CSOFT_CLIP = 4;
	/* ! @abstract CIGAR: H = clip on the read with clipped sequence trimmed off */
	public static final int BAM_CHARD_CLIP = 5;
	/* ! @abstract CIGAR: P = padding */
	public static final int BAM_CPAD = 6;
	/* ! @abstract CIGAR: equals = match */
	public static final int BAM_CEQUAL = 7;
	/* ! @abstract CIGAR: X = mismatch */
	public static final int BAM_CDIFF = 8;
	public static final int BAM_CBACK = 9;

	public static Map<Integer, Character> cigar2String = new HashMap<Integer, Character>(){
		{
			put(BAM_CMATCH, 'M');
			put(BAM_CINS, 'I');
			put(BAM_CDEL, 'D');
			put(BAM_CREF_SKIP, 'N');
			put(BAM_CSOFT_CLIP, 'S');
			put(BAM_CHARD_CLIP, 'H');
			put(BAM_CPAD, 'P');
			put(BAM_CEQUAL, '=');
			put(BAM_CDIFF, 'X');
		}
	};

	public static final int indelCallingWindowSize = 1000;

	/**
	 * 存储碱基信息（4bit）的个数
	 */
	private static final int CAPACITY;

	/**
	 * 二倍体基因型FASTA格式编码(A,C,G,T)
	 */
	private static final char[] FASTA_ABBR;

	/**
	 * VCF中基因型FASTA格式编码(A,C,G,T)
	 */
	private static final String[] FASTA_VABBR;

	/**
	 * 静态初始化块，用于加载属性文件
	 */
	static {
		CAPACITY = Byte.SIZE / 4;
		FASTA_ABBR = new char[] { 'A', 'M', 'W', 'R', 'M', 'C', 'Y', 'S', 'W',
				'Y', 'T', 'K', 'R', 'S', 'K', 'G', 'N' };
		FASTA_VABBR = new String[] { "A", "A,C", "A,T", "A,G", "A,C", "C",
				"C,T", "C,G", "A,T", "C,T", "T", "G,T", "A,G", "C,G", "G,T",
				"G", "N" };
	}

	/**
	 * 获取长整形变量的存储空间大小
	 * 
	 * @param
	 * @return int
	 */
	public static int getCapacity() {
		return CAPACITY;
	}

	/**
	 * 获取FASTA格式的字母表
	 * 
	 * @param
	 * @return char
	 */
	public static char getFastaAbbr(int index) {
		return FASTA_ABBR[index];
	}

	/**
	 * VCF获取FASTA格式的字母表
	 * 
	 * @param
	 * @return String
	 */
	public static String getFastaVabbr(int index) {
		return FASTA_VABBR[index];
	}

	public static final int VC_NO_GENO = 2;
	public static final int VC_BCFOUT = 4;
	public static final int VC_CALL = 8;
	public static final int VC_VARONLY = 16;
	public static final int VC_VCFIN = 32;
	public static final int VC_UNCOMP = 64;
	public static final int VC_KEEPALT = 256;
	public static final int VC_ACGT_ONLY = 512;
	public static final int VC_QCALL = 1024;
	public static final int VC_CALL_GT = 2048;
	public static final int VC_ADJLD = 4096;
	public static final int VC_NO_INDEL = 8192;
	public static final int VC_ANNO_MAX = 16384;
	public static final int VC_FIX_PL = 32768;
	public static final int VC_EM = 0x10000;
	public static final int VC_PAIRCALL = 0x20000;
	public static final int VC_QCNT = 0x40000;

	public static int B2B_INDEL_NULL = 10000;

	public static int RAND_MAX = 0x7fff;

	public static int DEF_MAPQ = 20;

	public static int CAP_DIST = 25;

	public static int BAM_CIGAR_SHIFT = 4;

	public static int BAM_CIGAR_MASK = (1 << BAM_CIGAR_SHIFT) - 1;

	public static int MINUS_CONST = 0x10000000;

	public static int INDEL_WINDOW_SIZE = 50;

	public static double M_LN10 = 2.30258509299404568402;

	public static double M_LN2 = 0.693147180559945309417;

	public static float CALL_DEFTHETA = 0.83f;

	public static int bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4,
			4, 4, 4, 4, 4 };

	public static int bam_nt16_table[] = { 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 1, 2, 4, 8, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			0 /* = */, 15, 15, 15, 1, 14, 2, 13, 15, 15, 4, 11, 15, 15, 12, 15,
			3, 15, 15, 15, 15, 5, 6, 8, 15, 7, 9, 15, 10, 15, 15, 15, 15, 15,
			15, 15, 1, 14, 2, 13, 15, 15, 4, 11, 15, 15, 12, 15, 3, 15, 15, 15,
			15, 5, 6, 8, 15, 7, 9, 15, 10, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15 };

	/**
	 * 单倍体基因型FASTA格式编码(A,C,G,T)
	 */
	public static final char[] FASTA_ABB;

	/**
	 * 静态初始化块，用于加载属性文件
	 */
	static {
		FASTA_ABB = new char[] { 'A', 'C', 'T', 'G' };
	}

	/**
	 * 获取单倍字母表
	 */
	public static char getFastaAbb(int code) {
		int index = (code & 0x07);
		if ((index & 0x04) != 0) {
			return 'N';
		} else
			return FASTA_ABB[index];
	}
}