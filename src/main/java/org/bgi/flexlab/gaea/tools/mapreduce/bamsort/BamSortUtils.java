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
package org.bgi.flexlab.gaea.tools.mapreduce.bamsort;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamFileHeaderMerger;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.seqdoop.hadoop_bam.SAMFormat;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class BamSortUtils {

    public static SAMFileHeader replaceSampleName(SAMFileHeader header,HashMap<String,String> replaceList){
        SAMFileHeader newHeader = header.clone();

        List<SAMReadGroupRecord> groups = header.getReadGroups();
        List<SAMReadGroupRecord> rgroups = new ArrayList<SAMReadGroupRecord>(
                groups.size());

        for (SAMReadGroupRecord g : groups) {
            if(!replaceList.containsKey(g.getSample()))
                throw new RuntimeException("replace list has no sample "+g.getSample());
            String newSample = replaceList.get(g.getSample());
            g.setSample(newSample);
            rgroups.add(g);
        }
        newHeader.setReadGroups(rgroups);
        return newHeader;
    }

    public static SAMFileHeader deleteSampleFromHeader(SAMFileHeader header, String sampleName) {
        SAMFileHeader newHeader = header.clone();
        List<SAMReadGroupRecord> newReadGroups = new ArrayList<>();

        for (SAMReadGroupRecord g : header.getReadGroups()) {
            if (g.getSample().equals(sampleName))
                newReadGroups.add(g);
        }
        newHeader.setReadGroups(newReadGroups);
        return newHeader;
    }

    public static String formatSampleName(String sampleName) {
        String formatSampleName;
        formatSampleName = sampleName.replaceAll("-", "");
        formatSampleName = formatSampleName.replaceAll("_", "");
        formatSampleName = formatSampleName.replaceAll("[*]", "x");
        return formatSampleName;
    }

    public static String getFileSuffix(SAMFormat format) {
        String fileSuffix = ".sorted.";
        if (format == null)
            fileSuffix += "cram";
        else if (format == SAMFormat.BAM)
            fileSuffix += "bam";
        else
            fileSuffix += "sam";
        return fileSuffix;
    }
}
