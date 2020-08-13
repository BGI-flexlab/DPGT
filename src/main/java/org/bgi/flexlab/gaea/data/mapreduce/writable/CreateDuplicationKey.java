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
package org.bgi.flexlab.gaea.data.mapreduce.writable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.RandomUtils;

public class CreateDuplicationKey {
	private String LB;
	private int chrIndex;
	private int position;
	private boolean forward;
	
	private final static String UNKNOW_LIBARY = "UNKNOW_LIBARY";
	private final SAMFileHeader header; 
	
	public CreateDuplicationKey(SAMFileHeader header) {
		this.header = header;
	}
	
	/**
	 * get duplication key;
	 */
	public void getKey(SAMRecord s,DuplicationKeyWritable key) {
		GaeaSamRecord sam = new GaeaSamRecord(header,s);
		sam.setHeader(header);
		sam.setDuplicateReadFlag(false);
		LB = getLibrary(sam);
					
		//unmaped
		if(sam.getReadUnmappedFlag()) {
			chrIndex = -1;
			position = RandomUtils.getRandomGenerator().nextInt()%100;
			key.set(LB, chrIndex, position, forward);
			return;
		}
			
		//SE && single-mapped
		if(!sam.getReadPairedFlag() || sam.getMateUnmappedFlag()) {
			chrIndex = sam.getReferenceIndex();
			position = sam.getReadNegativeStrandFlag() ? sam.getUnclippedEnd() : sam.getUnclippedStart();
			forward = !sam.getReadNegativeStrandFlag();
		}
		
		//both mapped
		if(sam.getReadPairedFlag() && !sam.getMateUnmappedFlag()) {
			if(sam.getReferenceIndex() > sam.getMateReferenceIndex() || (sam.getReferenceIndex() == sam.getMateReferenceIndex() && sam.getAlignmentStart() >= sam.getMateAlignmentStart())) {
				chrIndex = sam.getReferenceIndex();
				position = sam.getReadNegativeStrandFlag() ? sam.getUnclippedEnd() : sam.getUnclippedStart();
				forward = !sam.getReadNegativeStrandFlag();
			} else {
				chrIndex = sam.getMateReferenceIndex();
				position = (Integer) (sam.getMateNegativeStrandFlag() ? sam.getAttribute("ME") : sam.getAttribute("MS"));
				forward = !sam.getMateNegativeStrandFlag();
			}
		}
		
		key.set(LB, chrIndex, position, forward);
	}
	
	/**
	 * get library ID from RG of SAM record
	 */
	private String getLibrary(SAMRecord sam) {
		String Library;
		String readgroup = (String)sam.getAttribute("RG");
		if (readgroup == null) {
			Library = UNKNOW_LIBARY;		
		} else {
			Library=header.getReadGroup(readgroup).getLibrary();
			if (Library == null) {
				Library = header.getReadGroup(readgroup).getSample();
			}
		}
		return Library;
	}
}
