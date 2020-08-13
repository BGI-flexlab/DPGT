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
package org.bgi.flexlab.gaea.data.structure.header;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;

import java.util.ArrayList;
import java.util.List;

public class SamFileHeader {
	protected static boolean contains(SAMFileHeader header,
			ArrayList<SAMFileHeader> list) {
		for (SAMFileHeader that : list) {
			if (header.equals(that))
				return true;
		}
		return false;
	}
	
	public static SAMFileHeader createHeaderFromSampleName(
			SAMFileHeader header, String sampleName) {
		SAMFileHeader newHeader = header.clone();
		List<SAMReadGroupRecord> groups = header.getReadGroups();
		List<SAMReadGroupRecord> rgroups = new ArrayList<SAMReadGroupRecord>(
				groups.size());

		for (SAMReadGroupRecord g : groups) {
			if (g.getSample().equals(sampleName))
				rgroups.add(g);
		}
		newHeader.setReadGroups(rgroups);
		return newHeader;
	}
}
