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
package org.bgi.flexlab.gaea.data.structure.bam.filter;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.bam.filter.util.ReadsFilter;
import org.bgi.flexlab.gaea.util.ReadUtils;

import htsjdk.samtools.SAMRecord;

public class RealignerFilter extends ReadsFilter {

	public static boolean needToClean(GaeaSamRecord read,int maxInsertSize) {
		return !(read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag()
				|| read.getReadFailsVendorQualityCheckFlag()
				|| read.getMappingQuality() == 0
				|| read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START
				|| iSizeTooBigToMove(read, maxInsertSize)
				|| ReadUtils.is454Read(read) || ReadUtils.isIonRead(read));
	}

	public static boolean iSizeTooBigToMove(GaeaSamRecord read,
			int maxInsertSizeForMovingReadPairs) {
		return (read.getReadPairedFlag() && !read.getMateUnmappedFlag() && !read
				.getReferenceName().equals(read.getMateReferenceName())) 
				|| Math.abs(read.getInferredInsertSize()) > maxInsertSizeForMovingReadPairs; 
	}
}
