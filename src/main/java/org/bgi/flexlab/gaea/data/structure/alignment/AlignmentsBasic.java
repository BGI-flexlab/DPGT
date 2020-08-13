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
package org.bgi.flexlab.gaea.data.structure.alignment;

import htsjdk.samtools.SAMReadGroupRecord;
import org.bgi.flexlab.gaea.data.structure.bam.SAMCompressionInformationBasic;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AlignmentsBasic extends SAMCompressionInformationBasic implements Cloneable, Serializable {
	protected static Map<Integer, String> Id2Sample = new HashMap<Integer, String>();
	protected static Map<String, Integer> Sample2Id = new HashMap<String, Integer>();

	protected int sampleIndex;

	public AlignmentsBasic() {
		super();
	}

	public AlignmentsBasic( AlignmentsBasic read) {
		super(read);
		this.sampleIndex = read.sampleIndex;
	}

	public static void initIdSampleHash(List<SAMReadGroupRecord> samReadGroupRecords) {
		int i = 0;
		for(SAMReadGroupRecord samReadGroupRecord : samReadGroupRecords) {
			Id2Sample.put(i, samReadGroupRecord.getSample());
			Sample2Id.put(samReadGroupRecord.getSample(), i);
		}
	}

	public String getSample() {
		return Id2Sample.get(sampleIndex);
	}

	public int getSampleIndex() {
		return sampleIndex;
	}

	public void setSampleIndex(int sampleIndex) {
		this.sampleIndex = sampleIndex;
	}

	public void setSampleIndex(String sample) {
		this.sampleIndex = Sample2Id.get(sample);
	}
}

