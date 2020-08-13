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
package org.bgi.flexlab.gaea.framework.tools.mapreduce;

import htsjdk.samtools.SAMRecord;
import org.bgi.flexlab.gaea.data.mapreduce.writable.AlignmentBasicWritable;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;

/**
 * Created by zhangyong on 2017/3/8.
 */
public class WindowsBasedAlignmentMapper extends WindowsBasedMapper<AlignmentBasicWritable>{
    @Override
    void otherSetup(Context context) {
        AlignmentsBasic.initIdSampleHash(header.getReadGroups());
    }

    @Override
    void setOutputValue(SAMRecord samRecord) {
        AlignmentsBasic alignmentsBasic = new AlignmentsBasic();
        samRecord.setHeader(header);
        alignmentsBasic.parseSAM(samRecord);
        alignmentsBasic.setSampleIndex(samRecord.getReadGroup().getSample());
        outputValue.setAlignment(alignmentsBasic);
    }

    @Override
    void initOutputVaule() {
        outputValue = new AlignmentBasicWritable();
    }
}
