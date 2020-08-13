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

import org.bgi.flexlab.gaea.data.mapreduce.input.fastq.FastqMultipleSample;
import org.bgi.flexlab.gaea.data.structure.reads.ReadBasicStatistics;
import org.bgi.flexlab.gaea.data.structure.reads.ReadInformationWithSampleID;
import org.bgi.flexlab.gaea.data.structure.reads.report.FastqQualityControlReport;
import org.bgi.flexlab.gaea.tools.mapreduce.fastqqualitycontrol.FastqQualityControlOptions;

import java.io.IOException;
import java.util.ArrayList;

public class FastqQualityControlFilter {
	private FastqQualityControlOptions option;
	private AdaptorDynamicFilter adpFilter;
	private boolean dycut = false;
	private FastqQualityControlReport report;
	private int qualTrimNum = 0;

	public FastqQualityControlFilter(FastqQualityControlOptions option) {
		this.option = option;
		if (option.isQualityTrim())
			qualTrimNum = option.getQuality() - 33;
		adpFilter = new AdaptorDynamicFilter(option);

		int sampleSize = 1;
		if (option.getMultiSampleList() != null
				&& option.getMultiSampleList() != "") {
			FastqMultipleSample multiSampleList = null;
			try {
				multiSampleList = new FastqMultipleSample(
						option.getMultiSampleList(), false);
			} catch (IOException e) {
				e.printStackTrace();
			}
			sampleSize = multiSampleList.getSampleNumber();
		}
		report = new FastqQualityControlReport(sampleSize,
				option.isMultiStatis());
	}

	public boolean isDynamicCutted() {
		return dycut;
	}

	public ArrayList<ReadInformationWithSampleID> parseLine(
			ReadBasicStatistics stat, ArrayList<String> values) {
		ArrayList<ReadInformationWithSampleID> list = new ArrayList<ReadInformationWithSampleID>();

		for (String value : values) {
			String[] readLine = value.split("\t");
			if (readLine.length == 1) {
				stat.setAdaptor(Integer.parseInt(readLine[0]) - 1);
			} else {
				ReadInformationWithSampleID fqInfo = new ReadInformationWithSampleID(
						readLine);
				int[] trim = option.getTrim();
				if(fqInfo.isFirstRead()){
					if (trim[0] != 0 || trim[1] != 0) {
						fqInfo.trim(trim[0], trim[1]);
					}
				}else {
					if (trim[2] != 0 || trim[3] != 0) {
						fqInfo.trim(trim[2], trim[3]);
					}
				}
				int index = readLine[0].lastIndexOf("/");
				String tempkey = readLine[0].substring(index + 1).trim();
				stat.countBase(fqInfo, Integer.parseInt(tempkey) - 1);
				list.add(fqInfo);
			}
		}

		return list;
	}

	public boolean isBadReads(ReadBasicStatistics stat,
			ArrayList<ReadInformationWithSampleID> list) {
		if (stat.getReadCount() == 0) {
			return true;
		} else if (stat.getReadCount() == 1) {
			if (!option.isSEdata() && option.isFilterSE())
				return true;
		} else if (stat.getReadCount() > 2) {
			throw new RuntimeException(
					"more than 2 reads under same readsID without 1,2.\n"
							+ list.get(0).getReadName());
		}

		if (stat.isSampleIDException())
			throw new RuntimeException("pair read has difference sample id.\n"
					+ list.get(0).getReadName());

		if (list.size() > 1) {
			if (list.get(0).equals(list.get(1))) {
				throw new RuntimeException(
						"more than 2 reads under same readsID.\n"
								+ list.get(0).getReadName());
			}
		}
		return false;
	}

	public String qualityControlFilter(ArrayList<String> values) {
		ReadBasicStatistics stat = new ReadBasicStatistics(option);
		StringBuilder consensusSeq = new StringBuilder();
		ArrayList<ReadInformationWithSampleID> list = parseLine(stat, values);

		if (isBadReads(stat, list))
			return null;

		int sampleID = 0;
		if (!stat.getSampleID().equals("+"))
			sampleID = Integer.parseInt(stat.getSampleID());

		report.countRawReadInfo(stat, sampleID);

		boolean isClean = stat.getProblemReadsNum() == 0 ? true : false;
		report.countBaseByPosition(stat, sampleID, isClean);
		if (isClean) {
			report.countCleanReadInfo(stat, sampleID);
			int cnt = 0;
			for (ReadInformationWithSampleID read : list) {
				cnt++;

				consensusSeq.append(read.toString(qualTrimNum));
				if (cnt != list.size())
					consensusSeq.append("\n");
			}
		}else
			return null;

		return consensusSeq.toString();
	}

	public String filter(ArrayList<String> reads) {
		ArrayList<String> dyFilterReads = new ArrayList<String>();
		if (option.isDyncut()) {
			dycut = adpFilter.dyncutFilter(reads, dyFilterReads);
			if (dyFilterReads.size() == 0) {
				return null;
			}
		} else
			dyFilterReads = reads;

		return qualityControlFilter(dyFilterReads);
	}

	public FastqQualityControlReport getReport() {
		return this.report;
	}
}
