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
package org.bgi.flexlab.gaea.data.structure.bam.filter.util;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.Iterator;

/**
 * Created by zhangyong on 2017/3/9.
 */
public class FiltersMethod {

    public static boolean filterUnmappedReads(SAMRecord read) {
        return read.getReadUnmappedFlag() || read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
    }

    public static boolean filterNotPrimaryAlignment(SAMRecord read) {
        return read.getNotPrimaryAlignmentFlag();
    }

    public static boolean filterDuplicateRead(SAMRecord read) {
        return read.getDuplicateReadFlag();
    }

    public static boolean filterMappingQualityUnavailable(SAMRecord read,int unavailableQuality) {
        return (read.getMappingQuality() == unavailableQuality);
    }

    public static boolean FailsVendorQualityCheckFilter(SAMRecord read) {
        return read.getReadFailsVendorQualityCheckFlag();
    }

    public static boolean filterBadCigar(SAMRecord read) {
        final Cigar c = read.getCigar();

        // if there is no Cigar then it can't be bad
        if (c.isEmpty()) {
            return false;
        }

        Iterator<CigarElement> elementIterator = c.getCigarElements()
                .iterator();

        CigarOperator firstOp = CigarOperator.H;
        while (elementIterator.hasNext()
                && (firstOp == CigarOperator.H || firstOp == CigarOperator.S)) {
            CigarOperator op = elementIterator.next().getOperator();

            // No reads with Hard/Soft clips in the middle of the cigar
            if (firstOp != CigarOperator.H && op == CigarOperator.H) {
                return true;
            }
            firstOp = op;
        }

        // No reads starting with deletions (with or without preceding clips)
        if (firstOp == CigarOperator.D) {
            return true;
        }

        boolean hasMeaningfulElements = (firstOp != CigarOperator.H && firstOp != CigarOperator.S);
        boolean previousElementWasIndel = firstOp == CigarOperator.I;
        CigarOperator lastOp = firstOp;
        CigarOperator previousOp = firstOp;

        while (elementIterator.hasNext()) {
            CigarOperator op = elementIterator.next().getOperator();

            if (op != CigarOperator.S && op != CigarOperator.H) {

                // No reads with Hard/Soft clips in the middle of the cigar
                if (previousOp == CigarOperator.S
                        || previousOp == CigarOperator.H)
                    return true;

                lastOp = op;

                if (!hasMeaningfulElements && op.consumesReadBases()) {
                    hasMeaningfulElements = true;
                }

                if (op == CigarOperator.I || op == CigarOperator.D) {

                    // No reads that have consecutive indels in the cigar (II,
                    // DD, ID or DI)
                    if (previousElementWasIndel) {
                        return true;
                    }
                    previousElementWasIndel = true;
                } else {
                    previousElementWasIndel = false;
                }
            }
            // No reads with Hard/Soft clips in the middle of the cigar
            else if (op == CigarOperator.S && previousOp == CigarOperator.H) {
                return true;
            }

            previousOp = op;
        }

        // No reads ending in deletions (with or without follow-up clips)
        // No reads that are fully hard or soft clipped
        return lastOp == CigarOperator.D || !hasMeaningfulElements;
    }
    
    public static boolean matchBaseNumberZero(SAMRecord read){
    	
    	for(CigarElement ce : read.getCigar().getCigarElements()){
    		CigarOperator op = ce.getOperator();
    		if(op == CigarOperator.EQ || op == CigarOperator.M || op == CigarOperator.X){
    			return false;
    		}
    	}
    	
    	return true;
    }

    public static boolean filterBadMate(SAMRecord read) {
        return (read.getReadPairedFlag() && !read.getMateUnmappedFlag() && !read
                .getReferenceIndex().equals(read.getMateReferenceIndex()));
    }
}
