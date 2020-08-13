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
import org.bgi.flexlab.gaea.tools.mapreduce.fastqqualitycontrol.FastqQualityControlOptions;

import java.util.ArrayList;

public class AdaptorDynamicFilter {
	private DynamicAdaptor dpt;
	private int[] prep = null;
	private int adaptorLength = 10;

	public AdaptorDynamicFilter(FastqQualityControlOptions option) {
		adaptorLength = option.getAdaptorLength();

		dpt = new DynamicAdaptor();
		dpt.setArgs(option.getCest(), option.getBias(), adaptorLength,
				option.getMinimum());
		prep = dpt.BMprep(dpt.illumina, adaptorLength);
	}

	public boolean dyncutFilter(ArrayList<String> values, ArrayList<String> list) {
		ReadInformationWithSampleID read1 = new ReadInformationWithSampleID();
		ReadInformationWithSampleID read2 = new ReadInformationWithSampleID();
		boolean adp = false, r1 = false, r2 = false;
		for (String value : values) {
			String[] reads = value.split("\t");
			if (reads.length == 1) {
				adp = true;
				list.add(reads[0]);
			} else {
				if (reads[0].charAt(reads[0].length() - 1) == '1') {
					r1 = true;
					read1.set(reads);
				} else {
					r2 = true;
					read2.set(reads);
				}
			}
		}
		if (adp || !r1 || !r2) {
			if (r1)
				list.add(read1.getReadName() + "\t" + read1.getReadsSequence()
						+ "\t" + read1.getSampleID() + "\t"
						+ read1.getQualityString());
			if (r2)
				list.add(read2.getReadName() + "\t" + read2.getReadsSequence()
						+ "\t" + read2.getSampleID() + "\t"
						+ read2.getQualityString());
			return false;
		}
		int n = dpt
				.cut_adaptor(read1, read2, dpt.illumina, adaptorLength, prep);
		if ((n >> 31) != 0)
			return false;
		if (n != 0)
			return true;
		list.add(read1.getReadName() + "\t" + read1.getReadsSequence() + "\t"
				+ read1.getSampleID() + "\t" + read1.getQualityString());
		list.add(read2.getReadName() + "\t" + read2.getReadsSequence() + "\t"
				+ read2.getSampleID() + "\t" + read2.getQualityString());
		return false;
	}
}
