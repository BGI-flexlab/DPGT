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

import htsjdk.samtools.SAMFileHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.bgi.flexlab.gaea.data.mapreduce.input.bam.GaeaAnySAMInputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.partitioner.WindowsBasedComparator;
import org.bgi.flexlab.gaea.data.mapreduce.partitioner.WindowsBasedPartitioner;
import org.bgi.flexlab.gaea.data.mapreduce.partitioner.WindowsBasedSort;
import org.bgi.flexlab.gaea.data.structure.bam.filter.util.SamRecordFilter;
import org.seqdoop.hadoop_bam.SAMFormat;

import java.io.IOException;
import java.util.List;

public class
BioJob extends Job {

	@SuppressWarnings("deprecation")
	private BioJob(Configuration conf) throws IOException {
		super(conf);
	}

	public static BioJob getInstance() throws IOException {
		return getInstance(new Configuration());
	}

	public static BioJob getInstance(Configuration conf) throws IOException {
		return new BioJob(conf);
	}

	@SuppressWarnings("rawtypes")
	public void setWindowsBasicMapperClass(Class<? extends Mapper> cls, int windowsSize) throws IllegalStateException {
		conf.setInt(WindowsBasedMapper.WINDOWS_SIZE, windowsSize);
		setMapperClass(cls);
		setPartitionerClass(WindowsBasedPartitioner.class);
		setGroupingComparatorClass(WindowsBasedComparator.class);
		setSortComparatorClass(WindowsBasedSort.class);
		// setSortComparatorClass(WindowsBasedBasicSort.class);
	}

	@SuppressWarnings("rawtypes")
	public void setWindowsBasicMapperClass(Class<? extends Mapper> cls, int windowsSize, int windowsExtendSize)
			throws IllegalStateException {
		conf.setInt(WindowsBasedMapper.WINDOWS_EXTEND_SIZE, windowsExtendSize);
		setWindowsBasicMapperClass(cls, windowsSize);
	}

	public void setOnlyBaseRecalibrator(boolean bqsr) {
		conf.setBoolean(WindowsBasedMapper.BASERECALIBRATOR_ONLY, bqsr);
	}

	public void setWindowsBasicMapperClass(@SuppressWarnings("rawtypes") Class<? extends Mapper> cls, int windowsSize,
			int windowsExtendSize, boolean multiSample) throws IllegalStateException {
		if (multiSample)
			setMultipleSample();
		setWindowsBasicMapperClass(cls, windowsSize, windowsExtendSize);
	}

	public void setMultipleSample() {
		conf.setBoolean(WindowsBasedMapper.MULTIPLE_SAMPLE, true);
	}

	/*
	 * set windows based mapper filter
	 */
	public void setFilterClass(Class<? extends SamRecordFilter> cls) {
		conf.setClass(WindowsBasedMapper.SAM_RECORD_FILTER, cls, SamRecordFilter.class);
	}

	public void setOutputKeyValue(Class<?> mapKeyClass, Class<?> mapValueClass, Class<?> reduceKeyClass,
			Class<?> reduceValueClass) {
		setMapOutputKeyClass(mapKeyClass);
		setMapOutputValueClass(mapValueClass);
		setOutputKeyClass(reduceKeyClass);
		setOutputValueClass(reduceValueClass);
	}

	public void setOutputKeyValue(Class<?> keyClass, Class<?> valueClass) {
		setOutputKeyClass(keyClass);
		setOutputValueClass(valueClass);
	}

	/*
	 * set sam or bam inputformat
	 */
	public void setAnySamInputFormat(SAMFormat fmt) {
		conf.set(GaeaAnySAMInputFormat.SAM_FORMAT_FOR_ALL_PATH, fmt.toString());
		setInputFormatClass(GaeaAnySAMInputFormat.class);
	}
	
	public SAMFileHeader setHeader(List<Path> inputs , Path output) {
		try {
			return SamHdfsFileHeader.loadHeader(inputs, conf, output,false);
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}

	public SAMFileHeader setHeader(Path input, Path output) {
		try {
			return SamHdfsFileHeader.loadHeader(input, conf, output);
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}

	public void setHeader(String headerPath) {
		conf.set(SamHdfsFileHeader.BAM_HEADER_FILE_NAME, headerPath);
	}
}
