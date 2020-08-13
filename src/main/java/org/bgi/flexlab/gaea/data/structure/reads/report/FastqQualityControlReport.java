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
package org.bgi.flexlab.gaea.data.structure.reads.report;

import org.bgi.flexlab.gaea.data.structure.reads.ReadBasicStatistics;
import org.bgi.flexlab.gaea.util.ArrayListLongWrap;

import java.math.RoundingMode;
import java.text.DecimalFormat;

public class FastqQualityControlReport {
	public final static int STATIC_COUNT = 24;
	public final static int BASE_STATIC_COUNT = 14;
	public final static int BASE_START_INDEX = STATIC_COUNT - BASE_STATIC_COUNT;

	private long[][] statNumbers;

	private ArrayListLongWrap[][] BASE_BY_POSITION;

	private int sampleSize;
	private boolean isMulti = true;

	public FastqQualityControlReport(int ssize, boolean isMulti) {
		this.isMulti = isMulti;

		sampleSize = ssize;
		if (!isMulti)
			sampleSize = 1;

		statNumbers = new long[sampleSize][STATIC_COUNT];

		BASE_BY_POSITION = new ArrayListLongWrap[sampleSize][BASE_STATIC_COUNT];
		initArrayListLongWrap(BASE_BY_POSITION);
	}

	private void initArrayListLongWrap(ArrayListLongWrap[][] data) {
		for (int i = 0; i < sampleSize; i++) {
			for (int j = 0; j < BASE_STATIC_COUNT; j++) {
				data[i][j] = new ArrayListLongWrap();
			}
		}
	}

	public void countRawReadInfo(ReadBasicStatistics stat, int sampleID) {
		if (!isMulti) {
			sampleID = 0;
		}
		statNumbers[sampleID][Type.TOTAL_READS_NUM.get()] += stat
				.getReadCount();
		statNumbers[sampleID][Type.RAW_TOTAL_BASE_NUM.get()] += stat
				.getBaseNumber();

		for (int i = 0; i < (BASE_STATIC_COUNT / 2); i++) {
			statNumbers[sampleID][i + BASE_START_INDEX] += stat
					.getBasicBaseCount(i);
		}

		if (stat.getProblemReadsNum() != 0) {
			statNumbers[sampleID][Type.ADAPTER_READS_NUM.get()] += stat
					.getAdaptorCount();
			statNumbers[sampleID][Type.TOO_MANY_N_READS_NUM.get()] += stat
					.getTooManyNCounter();
			statNumbers[sampleID][Type.TOO_MANY_LOW_QUAL_READS_NUM.get()] += stat
					.getLowQualityCounter();
			if (stat.getReadCount() == 2 && stat.getProblemReadsNum() == 1) {
				if (stat.getAdaptorCount() == 1)
					statNumbers[sampleID][Type.Adapter_FILTER_DUE_TO_PE.get()] += 1;
				else if (stat.getTooManyNCounter() == 1)
					statNumbers[sampleID][Type.tooManyN_FILTER_DUE_TO_PE.get()] += 1;
				else if (stat.getLowQualityCounter() == 1)
					statNumbers[sampleID][Type.lowQual_FILTER_DUE_TO_PE.get()] += 1;
			}
		}
	}

	public void countCleanReadInfo(ReadBasicStatistics stat, int sampleID) {
		if (!isMulti) {
			sampleID = 0;
		}
		statNumbers[sampleID][Type.CLEAN_TOTAL_READS_NUM.get()] += stat
				.getReadCount();
		statNumbers[sampleID][Type.CLEAN_TOTAL_BASE_NUM.get()] += stat
				.getBaseNumber();

		for (int i = (BASE_STATIC_COUNT / 2); i < BASE_STATIC_COUNT; i++) {
			statNumbers[sampleID][i + BASE_START_INDEX] += stat
					.getBasicBaseCount(i - (BASE_STATIC_COUNT / 2));
		}
	}

	public void countBaseByPosition(ReadBasicStatistics stat, int sampleID,
			boolean isClean) {
		if (!isMulti) {
			sampleID = 0;
		}

		int i;
		for (i = 0; i < (BASE_STATIC_COUNT / 2); i++) {
			BASE_BY_POSITION[sampleID][i].add(stat.getPositionInfo(i));
		}

		if (isClean) {
			int start = BASE_STATIC_COUNT / 2;
			for (i = start; i < BASE_STATIC_COUNT; i++) {
				BASE_BY_POSITION[sampleID][i].add(stat.getPositionInfo(i
						- start));
			}
		}
	}

	public long getCount(int sampleID, int index) {
		if (!isMulti)
			sampleID = 0;
		return statNumbers[sampleID][index];
	}

	public ArrayListLongWrap getBaseByPosition(int sampleID, int index) {
		if (!isMulti)
			sampleID = 0;
		return BASE_BY_POSITION[sampleID][index];
	}

	private boolean partitionNull = false;

	public boolean isPartitionNull() {
		return partitionNull;
	}

	public int addCount(String line) {
		String[] str = line.split("\t");
		int sampleID = Integer.parseInt(str[0]);

		if (!isMulti)
			sampleID = 0;

		for (int i = 1; i < str.length; i++) {
			statNumbers[sampleID][i - 1] += Long.parseLong(str[i]);
		}

		partitionNull = false;
		if (Long.parseLong(str[1]) == 0)
			partitionNull = true;

		return sampleID;
	}

	public void addBaseByPosition(int sampleID, int index, String line) {
		String[] str = line.split("\t");
		for (int i = 0; i < str.length; i++) {
			if (str[i].equals(""))
				continue;
			BASE_BY_POSITION[sampleID][index].add(i, Long.parseLong(str[i]));
		}
	}

	public String toString() {
		StringBuilder strBuilder = new StringBuilder();

		int i, j;
		for (i = 0; i < sampleSize; i++) {
			strBuilder.append(i);

			for (j = 0; j < FastqQualityControlReport.STATIC_COUNT; j++) {
				strBuilder.append("\t" + getCount(i, j));
			}

			if (getCount(i, 0) == 0) {
				continue;
			} else
				strBuilder.append("\n");

			for (j = 0; j < FastqQualityControlReport.BASE_STATIC_COUNT; j++) {
				strBuilder.append(getBaseByPosition(i, j).toString());
				if (j != (FastqQualityControlReport.BASE_STATIC_COUNT - 1))
					strBuilder.append("\n");
			}
		}
		return strBuilder.toString();
	}

	public enum Type {
		TOTAL_READS_NUM(0), CLEAN_TOTAL_READS_NUM(1), ADAPTER_READS_NUM(2), TOO_MANY_N_READS_NUM(
				3), TOO_MANY_LOW_QUAL_READS_NUM(4), Adapter_FILTER_DUE_TO_PE(5), tooManyN_FILTER_DUE_TO_PE(
				6), lowQual_FILTER_DUE_TO_PE(7), RAW_TOTAL_BASE_NUM(8), CLEAN_TOTAL_BASE_NUM(
				9), RAW_A_BASE_NUM(10), RAW_C_BASE_NUM(11), RAW_T_BASE_NUM(12), RAW_G_BASE_NUM(
				13), RAW_N_BASE_NUM(14), RAW_Q20_NUM(15), RAW_Q30_NUM(16), CLEAN_A_BASE_NUM(
				17), CLEAN_C_BASE_NUM(18), CLEAN_T_BASE_NUM(19), CLEAN_G_BASE_NUM(
				20), CLEAN_N_BASE_NUM(21), CLEAN_Q20_NUM(22), CLEAN_Q30_NUM(23);

		private int index;

		private Type(int i) {
			index = i;
		}

		public int get() {
			return index;
		}
	}

	public String getReportContext(int sampleID) {
		DecimalFormat df = new DecimalFormat("0.000");
		df.setRoundingMode(RoundingMode.HALF_UP);

		StringBuilder outString = new StringBuilder();
		outString.append("Filter Information:\n");
		outString.append("Total reads number: ");
		outString.append(statNumbers[sampleID][Type.TOTAL_READS_NUM.get()]);
		outString.append("\nClean reads number:");
		outString
				.append(statNumbers[sampleID][Type.CLEAN_TOTAL_READS_NUM.get()]);
		outString.append("\nAdapter reads number: ");
		outString.append(statNumbers[sampleID][Type.ADAPTER_READS_NUM.get()]);
		outString.append("\nToo many low quality reads number: ");
		outString.append(statNumbers[sampleID][Type.TOO_MANY_LOW_QUAL_READS_NUM
				.get()]);
		outString.append("\nToo many N bases reads number: ");
		outString
				.append(statNumbers[sampleID][Type.TOO_MANY_N_READS_NUM.get()]);
		outString
				.append("\nFiltered reads due to another PE reads is adapter: ");
		outString.append(statNumbers[sampleID][Type.Adapter_FILTER_DUE_TO_PE
				.get()]);
		outString
				.append("\nFiltered reads due to another PE reads has too many low quality bases: ");
		outString.append(statNumbers[sampleID][Type.lowQual_FILTER_DUE_TO_PE
				.get()]);
		outString
				.append("\nFiltered reads due to another PE reads has too many N bases: ");
		outString.append(statNumbers[sampleID][Type.tooManyN_FILTER_DUE_TO_PE
				.get()]);
		outString.append("\nFiltered reads persentage: ");
		outString
				.append(df.format(100
						* (statNumbers[sampleID][Type.TOTAL_READS_NUM.get()] - statNumbers[sampleID][Type.CLEAN_TOTAL_READS_NUM
								.get()])
						/ (double) statNumbers[sampleID][Type.TOTAL_READS_NUM
								.get()]));
		outString.append("%");
		outString.append("\n\nRaw Reads Information:\n");
		outString.append("Total base Number: ");
		outString.append(statNumbers[sampleID][Type.RAW_TOTAL_BASE_NUM.get()]);
		outString.append("\nGC rate: ");
		outString
				.append(df.format(100
						* (statNumbers[sampleID][Type.RAW_C_BASE_NUM.get()] + statNumbers[sampleID][Type.RAW_G_BASE_NUM
								.get()])
						/ (double) (statNumbers[sampleID][Type.RAW_TOTAL_BASE_NUM
								.get()] - statNumbers[sampleID][Type.RAW_N_BASE_NUM
								.get()])));
		outString.append("%\nN base rate: ");
		outString
				.append(df.format(100
						* statNumbers[sampleID][Type.RAW_N_BASE_NUM.get()]
						/ (double) statNumbers[sampleID][Type.RAW_TOTAL_BASE_NUM
								.get()]));
		outString.append("%\nQuality>=20 base rate: ");
		outString
				.append(df.format(100
						* statNumbers[sampleID][Type.RAW_Q20_NUM.get()]
						/ (double) statNumbers[sampleID][Type.RAW_TOTAL_BASE_NUM
								.get()]));
		outString.append("%\nQuality>=30 base rate: ");
		outString
				.append(df.format(100
						* statNumbers[sampleID][Type.RAW_Q30_NUM.get()]
						/ (double) statNumbers[sampleID][Type.RAW_TOTAL_BASE_NUM
								.get()]));
		outString.append("%\n\nClean Reads Information:\n");
		outString.append("Total base Number: ");
		outString
				.append(statNumbers[sampleID][Type.CLEAN_TOTAL_BASE_NUM.get()]);
		outString.append("\nGC rate: ");
		outString
				.append(df.format(100
						* (statNumbers[sampleID][Type.CLEAN_C_BASE_NUM.get()] + statNumbers[sampleID][Type.CLEAN_G_BASE_NUM
								.get()])
						/ (double) (statNumbers[sampleID][Type.CLEAN_TOTAL_BASE_NUM
								.get()] - statNumbers[sampleID][Type.CLEAN_N_BASE_NUM
								.get()])));
		outString.append("%\nN base rate: ");
		outString.append(df.format(100
				* statNumbers[sampleID][Type.CLEAN_N_BASE_NUM.get()]
				/ (double) statNumbers[sampleID][Type.CLEAN_TOTAL_BASE_NUM
						.get()]));
		outString.append("%\nQuality>=20 base rate: ");
		outString.append(df.format(100
				* statNumbers[sampleID][Type.CLEAN_Q20_NUM.get()]
				/ (double) statNumbers[sampleID][Type.CLEAN_TOTAL_BASE_NUM
						.get()]));
		outString.append("%\nQuality>=30 base rate: ");
		outString.append(df.format(100
				* statNumbers[sampleID][Type.CLEAN_Q30_NUM.get()]
				/ (double) statNumbers[sampleID][Type.CLEAN_TOTAL_BASE_NUM
						.get()]));
		outString.append("%\n\n\n");
		return outString.toString();
	}

	public String getGraphContext(int sampleID) {
		StringBuilder sb = new StringBuilder();
		int i;
		for (i = 0; i < 5; i++) {
			sb.append(BASE_BY_POSITION[sampleID][i].toString());
			sb.append("\n");
		}
		for (i = 7; i < 12; i++) {
			sb.append(BASE_BY_POSITION[sampleID][i].toString());
			sb.append("\n");
		}
		for (i = 5; i < 7; i++) {
			sb.append(BASE_BY_POSITION[sampleID][i].toString());
			sb.append("\n");
		}
		for (i = 12; i < 14; i++) {
			sb.append(BASE_BY_POSITION[sampleID][i].toString());
			sb.append("\n");
		}

		return sb.toString();
	}
}
