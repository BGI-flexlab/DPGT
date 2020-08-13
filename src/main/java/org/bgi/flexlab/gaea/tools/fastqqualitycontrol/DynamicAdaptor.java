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
package org.bgi.flexlab.gaea.tools.fastqqualitycontrol;

import org.bgi.flexlab.gaea.data.structure.reads.ReadInformationWithSampleID;

public class DynamicAdaptor {
	private int cest = 3;
	private int minimum = 0;
	private long filter_reads = 0;
	private int length_adaptor = 10;

	public final int[] illumina = { 2, 8, 4, 8, 2, 8, 2, 8, 8, 1, 8, 1, 2, 1,
			2, 1, 8, 2, 8 };
	private static final int[] dna_compTable = { 0, 8, 4, 12, 2, 10, 9, 14, 1,
			6, 5, 13, 3, 11, 7, 15 };
	private static final char dna_nt16Table[] = { 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 15, 15, 15, 15, 15, 1, 2, 4, 8, 15, 15, 15, 15, 15, 15, 15,
			15, 15, 0, 15, 15, 15, 1, 14, 2, 13, 15, 15, 4, 11, 15, 15, 12, 15,
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

	private static final int min32(int x, int y) {
		int min = y + (x - y & x - y >> 31);
		return min;
	}

	private String reverseSEQ(char[] seq, int len) {
		for (int i = 0; i < len >> 1; i++) {
			char tmp = seq[i];
			seq[i] = seq[(len - i - 1)];
			seq[(len - i - 1)] = tmp;
		}
		return new String(seq);
	}

	private void transSEQ(int[] s, int l) {
		int i;
		for (i = 0; i < l >> 1; i++) {
			int t = dna_compTable[s[(l - i - 1)]];
			s[(l - i - 1)] = dna_compTable[s[i]];
			s[i] = t;
		}
		if ((l & 0x1) != 0)
			s[i] = dna_compTable[s[i]];
	}

	private int[] seq2code(char[] seq, int length) {
		int[] c = new int[length];
		for (int i = 0; i < length; i++)
			c[i] = dna_nt16Table[seq[i]];
		return c;
	}

	public int[] BMprep(int[] pat, int m) {
		int i;
		int[] suff = (int[]) null;
		int[] prep = (int[]) null;
		prep = new int[m + 15];
		for (i = m; i < m + 15; i++)
			prep[i] = m;
		for (i = 0; i < m - 1; i++)
			prep[(pat[i] + m)] = (m - i - 1);
		suff = new int[m];
		int f = 0;
		int g = m - 1;
		for (i = m - 2; i >= 0; i--) {
			if ((i > g) && (suff[(i + m - 1 - f)] < i - g)) {
				suff[i] = suff[(i + m - 1 - f)];
			} else {
				if (i < g)
					g = i;
				f = i;
				while ((g >= 0) && (pat[g] == pat[(g + m - 1 - f)]))
					g--;
				suff[i] = (f - g);
			}
		}
		int j = 0;
		for (i = 0; i < m; i++)
			prep[i] = m;
		for (i = m - 1; i >= 0; i--)
			if (suff[i] == i + 1)
				for (; j < m - 1 - i; j++)
					prep[j] = (m - 1 - i);
		for (i = 0; i <= m - 2; i++)
			prep[(m - 1 - suff[i])] = (m - 1 - i);
		return prep;
	}

	private int location(int[] str, int n, int[] pat, int m, int[] prep) {
		int i, j = 0;
		while (j <= n - m) {
			int b = 0;
			int k = 0;
			for (i = m - 1; i >= 0; i--) {
				if ((pat[i] & str[(i + j)]) == 0) {
					if (b != 0) {
						i = k;
						break;
					}

					k = i;
					b++;
				}
			}

			if (i >= 0) {
				int max = prep[(str[(i + j)] + m)] - m + 1 + i;
				if (max < prep[i])
					max = prep[i];
				j += max;
			} else {
				return j;
			}
		}
		j = n - m;
		while (j < n - this.cest) {
			for (i = 0; (i < n - j) && (pat[i] == str[(j + i)]); i++)
				;
			if (i < n - j)
				j++;
			else
				return j;
		}
		return 0;
	}

	private void updateRead(ReadInformationWithSampleID r, int loc) {
		int len = r.getReadLength();
		String seq = reverseSEQ(r.getReadsSequence().toCharArray(), len);
		seq = seq.substring(0, loc);
		seq = reverseSEQ(seq.toCharArray(), loc);
		r.setReadSequence(seq);
		seq = reverseSEQ(r.getQualityString().toCharArray(), len);
		seq = seq.substring(0, loc);
		seq = reverseSEQ(seq.toCharArray(), loc);
		r.setReadQuality(seq);
	}

	public int cut_adaptor(ReadInformationWithSampleID r1,
			ReadInformationWithSampleID r2, int[] pat, int len, int[] prep) {
		int m = 0;
		int n = 0;
		int loc = 0;

		int[] s = seq2code(r1.getReadsSequence().toCharArray(),
				r1.getReadLength());
		int[] p = seq2code(r2.getReadsSequence().toCharArray(),
				r2.getReadLength());
		m = location(s, s.length, pat, len, prep);
		n = location(p, p.length, pat, len, prep);
		if ((m != 0) && (n != 0))
			loc = min32(m, n);
		if (loc != 0) {
			if ((this.minimum != 0) && (loc < this.minimum)) {
				this.filter_reads += 1L;
				return -1;
			}
			r1.setReadSequence(r1.getReadsSequence().substring(0, loc));
			r1.setReadQuality(r1.getQualityString().substring(0, loc));
			r2.setReadSequence(r2.getReadsSequence().substring(0, loc));
			r2.setReadQuality(r2.getQualityString().substring(0, loc));
		} else {
			transSEQ(s, r1.getReadLength());
			transSEQ(p, r2.getReadLength());
			m = location(s, s.length, pat, len, prep);
			n = location(p, p.length, pat, len, prep);
			if ((m != 0) && (n != 0))
				loc = min32(m, n);
			if (loc != 0) {
				if ((this.minimum != 0) && (loc < this.minimum)) {
					this.filter_reads += 1L;
					return -1;
				}
				updateRead(r1, loc);
				updateRead(r2, loc);
			} else {
				return 0;
			}
		}
		return 1;
	}

	public void setArgs(int t, int i, int s, int m) {
		this.cest = t;
		this.minimum = m;
		this.length_adaptor = s;
	}

	public long getFilterRead() {
		return this.filter_reads;
	}

	public int getLengthAdaptor() {
		return this.length_adaptor;
	}
}
