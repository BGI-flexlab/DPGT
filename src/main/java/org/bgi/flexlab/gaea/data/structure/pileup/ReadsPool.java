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
package org.bgi.flexlab.gaea.data.structure.pileup;

import htsjdk.samtools.SAMFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.AlignmentBasicWritable;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;

import java.util.Iterator;
/**
 * Created by zhangyong on 2016/12/23.
 */
public class ReadsPool implements Iterator {
    private boolean isSAM;

    private Iterator<GaeaSamRecord> samReads;

    private Iterator<AlignmentBasicWritable> alignments;

    public ReadsPool(Iterator<GaeaSamRecord> samReads, SAMFileHeader header) {
        this.samReads = samReads;
        isSAM = true;
    }

    public ReadsPool(Iterator<AlignmentBasicWritable> alignments) {
        this.alignments = alignments;
        isSAM = false;
    }

    @Override
    public AlignmentsBasic next() {
        AlignmentsBasic alignment;
        if(isSAM) {
            alignment = new AlignmentsBasic();
            alignment.parseSAM(samReads.next());
        } else {
            alignment = alignments.next().getAlignment();
        }

        return alignment;
    }

    @Override
    public boolean hasNext() {
        if (isSAM) {
            return samReads.hasNext();
        } else {
            return alignments.hasNext();
        }
    }
}
