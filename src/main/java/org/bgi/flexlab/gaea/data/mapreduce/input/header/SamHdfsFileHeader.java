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
package org.bgi.flexlab.gaea.data.mapreduce.input.header;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.structure.CramHeader;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.PathFilter;
import org.apache.hadoop.fs.permission.FsAction;
import org.apache.hadoop.fs.permission.FsPermission;
import org.bgi.flexlab.gaea.data.exception.FileNotExistException;
import org.bgi.flexlab.gaea.data.structure.header.SamFileHeader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/* *
 * bam header io for hdfs
 * */
public class SamHdfsFileHeader extends SamFileHeader {
	public final static String BAM_HEADER_FILE_NAME = "SAMFileHeader";
	protected final static SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;
	protected static boolean MERGE_SEQUENCE_DICTIONARIES = true;
	protected static SAMFileHeader.SortOrder headerMergerSortOrder;
	protected final static boolean ASSUME_SORTED = false;

	public static class HeaderPathFilter implements PathFilter {
		@Override
		public boolean accept(Path path) {
			if (path.getName().startsWith("_"))
				return false;
			return true;
		}
	}

	@SuppressWarnings("deprecation")
	public static SAMFileHeader getSAMHeader(FileSystem fs, Path file) throws IOException {
		SamReader samr = SamReaderFactory.makeDefault().open((java.nio.file.Path) file);
		SAMFileHeader header = samr.getFileHeader();
		samr.close();
		if (header == null) {
			throw new FileNotExistException.MissingHeaderException(file.getName());
		}

		return header;
	}

	public static CramHeader getCramHeader(FileSystem fs, Path path) throws IOException {
		SeekableStream sin = WrapSeekable.openPath(fs, path);
		CramHeader header = CramIO.readCramHeader(sin);
		sin.close();
		return header;
	}
	
	public static SAMFileHeader getSAMFileHeader(Path input, FileSystem fs,boolean cram) {
		SAMFileHeader header = null;
		if (!cram) {
			try {
				header = getSAMHeader(fs, input);
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}
		} else {
//			try {
//				header = getCramHeader(fs, input);
//			} catch (IOException e) {
//				throw new RuntimeException(e.toString());
//			}
		}
		return header;
	}
	
	/*public static SAMFileHeader traversal(Path input, FileSystem fs, Configuration conf, boolean cram) {
		FileStatus status = null;
		try {
			status = fs.getFileStatus(input);
		} catch (IOException e2) {
			throw new FileNotExistException(input.getName());
		}

		if (status.isFile()) {
			return getSAMFileHeader(input,fs,cram);
		}else {
			FileStatus[] stats = null;
			try {
				stats = fs.listStatus(input, new HeaderPathFilter());
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}
			List<Path> inputs = new ArrayList<Path>();
			for (FileStatus file : stats) {
				inputs.add(file.getPath());
			}
			
			return traversal(inputs,fs,conf,cram);
		}
	}*/

	public static SAMFileHeader traversal(List<Path> inputs , FileSystem fs ,Configuration conf , boolean cram) {
		ArrayList<SAMFileHeader> mergeHeaders = new ArrayList<SAMFileHeader>();
		SAMFileHeader mergedHeader = null;
		boolean matchedSortOrders = true;
		
		for(Path path : inputs) {
			SAMFileHeader header = getSAMFileHeader(path,fs,cram);
			matchedSortOrders = matchedSortOrders && header.getSortOrder() == SORT_ORDER;
			if (!contains(header, mergeHeaders))
				mergeHeaders.add(header);
		}
		
		if (matchedSortOrders || SORT_ORDER == SAMFileHeader.SortOrder.unsorted || ASSUME_SORTED) {
			headerMergerSortOrder = SORT_ORDER;
		} else {
			headerMergerSortOrder = SAMFileHeader.SortOrder.unsorted;
		}
		mergedHeader = new SamFileHeaderMerger(headerMergerSortOrder, mergeHeaders, MERGE_SEQUENCE_DICTIONARIES)
				.getMergedHeader();
		return mergedHeader;
	}

	public static SAMFileHeader traversal(Path input, FileSystem fs, Configuration conf, boolean cram) {
		ArrayList<SAMFileHeader> mergeHeaders = new ArrayList<SAMFileHeader>();
		SAMFileHeader mergedHeader = null;
		boolean matchedSortOrders = true;

		FileStatus status = null;
		try {
			status = fs.getFileStatus(input);
		} catch (IOException e2) {
			throw new FileNotExistException(input.getName());
		}

		if (status.isFile()) {
			SAMFileHeader header = getSAMFileHeader(input,fs,cram);
			matchedSortOrders = matchedSortOrders && header.getSortOrder() == SORT_ORDER;
			if (!contains(header, mergeHeaders))
				mergeHeaders.add(header);
		} else {
			FileStatus[] stats = null;
			try {
				stats = fs.listStatus(input, new HeaderPathFilter());
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}

			for (FileStatus file : stats) {
				Path filePath = file.getPath();
				SAMFileHeader header = null;
				if (file.isFile()) {
					header = getSAMFileHeader(filePath,fs,cram);
				} else {
					header = traversal(filePath, fs, conf, cram);
				}
				matchedSortOrders = matchedSortOrders && header.getSortOrder() == SORT_ORDER;
				if (!contains(header, mergeHeaders))
					mergeHeaders.add(header);
			}
		}
		if (matchedSortOrders || SORT_ORDER == SAMFileHeader.SortOrder.unsorted || ASSUME_SORTED) {
			headerMergerSortOrder = SORT_ORDER;
		} else {
			headerMergerSortOrder = SAMFileHeader.SortOrder.unsorted;
		}
		mergedHeader = new SamFileHeaderMerger(headerMergerSortOrder, mergeHeaders, MERGE_SEQUENCE_DICTIONARIES)
				.getMergedHeader();
		return mergedHeader;
	}

	public static SAMFileHeader loadHeader(Path input, Configuration conf, Path output) throws IOException {
		return loadHeader(input, conf, output, false);
	}

	public static SAMFileHeader loadHeader(Path input, Configuration conf, Path output, boolean cram)
			throws IOException {
//		List<Path> inputs = new ArrayList<Path>();
//		inputs.add(input);
//		return loadHeader(inputs,conf,output,cram);
		SAMFileHeader mergeHeader = traversal(input, input.getFileSystem(conf), conf, cram);
		writeHeader(conf, mergeHeader, output);
		return mergeHeader;
	}
	
	public static SAMFileHeader loadHeader(List<Path> inputs, Configuration conf, Path output, boolean cram)
			throws IOException {
		FileSystem fs = inputs.get(0).getFileSystem(conf);
		SAMFileHeader mergeHeader = traversal(inputs, fs, conf, cram);
		if (mergeHeader == null) {
			throw new FileNotExistException.MissingHeaderException(inputs.get(0).getName());
		}
		writeHeader(conf, mergeHeader, output);

		return mergeHeader;
	}

	public static void writeHeader(Configuration conf, SAMFileHeader header, Path output) {
		Path rankSumTestObjPath = null;
		FsAction[] v = FsAction.values();
		StringBuilder uri = new StringBuilder();
		uri.append(output);
		if (!output.getName().endsWith("/")) {
			uri.append(System.getProperty("file.separator"));
		}
		uri.append(BAM_HEADER_FILE_NAME);
		conf.set(BAM_HEADER_FILE_NAME, uri.toString());
		rankSumTestObjPath = new Path(uri.toString());
		FileSystem fs = null;
		try {
			fs = rankSumTestObjPath.getFileSystem(conf);
			FsPermission permission = new FsPermission(v[7], v[7], v[7]);
			if (!fs.exists(output)) {
				fs.mkdirs(output, permission);
			} else {
				fs.setPermission(output, permission);
			}

			SamFileHeaderCodec.writeHeader(header, fs.create(rankSumTestObjPath));
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		} finally {
			try {
				fs.close();
			} catch (IOException ioe) {
				throw new RuntimeException(ioe.toString());
			}
		}
	}

	public static SAMFileHeader getHeader(Configuration conf) {
		if (conf.get(BAM_HEADER_FILE_NAME) == null)
			return null;

		SAMFileHeader header = null;
		try {
			Path headerPath = new Path(conf.get(BAM_HEADER_FILE_NAME));
			HdfsHeaderLineReader reader = new HdfsHeaderLineReader(headerPath, conf);
			header = SamFileHeaderCodec.readHeader(reader);
		} catch (IOException e) {
			throw new RuntimeException(e.getMessage());
		}

		return header;
	}
}
