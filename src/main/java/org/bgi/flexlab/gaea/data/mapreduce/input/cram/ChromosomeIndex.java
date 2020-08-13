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
package org.bgi.flexlab.gaea.data.mapreduce.input.cram;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.cram.build.CramContainerIterator;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.LineReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class ChromosomeIndex {
	private String cramFileName;
	private SAMFileHeader samFileHeader;
	private String cramIndexFileName = null;

	private HashMap<Integer, ChromosomeIndexStructure> chromosomeIndexMap = null;

	class ChromosomeIndexStructure {
		private int sequenceId;
		private long start;
		private long end;

		public ChromosomeIndexStructure(int sid, long s, long e) {
			sequenceId = sid;
			start = s;
			end = e;
		}

		public ChromosomeIndexStructure(String line) {
			String[] str = line.split("\t");
			sequenceId = Integer.parseInt(str[0]);
			start = Long.parseLong(str[1]);
			end = Long.parseLong(str[2]);
		}

		public int getSeqId() {
			return sequenceId;
		}

		public long getStart() {
			return this.start;
		}

		public long getEnd() {
			return this.end;
		}

		public String toString() {
			return sequenceId + "\t" + start + "\t" + end;
		}
	}

	public ChromosomeIndex(String cram) {
		this.cramFileName = cram;
		cramIndexFileName = this.cramFileName + ".crai";
	}

	public ChromosomeIndex(String cram_, String index_) {
		this.cramFileName = cram_;
		this.cramIndexFileName = index_;
	}

	public void setHeader(SAMFileHeader header) {
		samFileHeader = header;
	}

	public ArrayList<ChromosomeIndexStructure> indexForChromosome(Path p) {
		long offset = 0;
		Container c = null;
		int prevSeqID = Integer.MIN_VALUE;
		long start = -1, end = -1;

		ArrayList<ChromosomeIndexStructure> chrIndex = new ArrayList<ChromosomeIndexStructure>();

		SeekableStream seekableStream = null;
		try {
			seekableStream = WrapSeekable.openPath(
					new Configuration(), p);
		} catch (IOException e) {
			e.printStackTrace();
		}
		CramContainerIterator cci = null;
		try {
			cci = new CramContainerIterator(seekableStream);
		} catch (IOException e) {
			e.printStackTrace();
		}

		samFileHeader = cci.getCramHeader().getSamFileHeader();

		try {
			offset = seekableStream.position();
		} catch (IOException e) {
			e.printStackTrace();
		}

		while (cci.hasNext()) {
			c = cci.next();
			c.offset = offset;
			if (c.sequenceId != prevSeqID) {
				if (prevSeqID != Integer.MIN_VALUE) {
					end = c.offset - 1;
					chrIndex.add(new ChromosomeIndexStructure(prevSeqID,
							start, end));
				}
				start = offset;
				prevSeqID = c.sequenceId;
			}
			try {
				offset = seekableStream.position();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		chrIndex.add(new ChromosomeIndexStructure(prevSeqID, start, end));

		return chrIndex;
	}

	public void writeIndex(ArrayList<ChromosomeIndexStructure> chrIndex) {
		Path p = new Path(cramIndexFileName);
		Configuration conf = new Configuration();
		FSDataOutputStream output = null;
		try {
			output = p.getFileSystem(conf).create(p);
		} catch (IOException e) {
			e.printStackTrace();
		}

		for (ChromosomeIndexStructure chrindex : chrIndex) {
			try {
				output.write(chrindex.toString().getBytes());
				output.write("\n".getBytes());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		try {
			output.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void loadChromosomeIndex() {
		chromosomeIndexMap = new HashMap<Integer, ChromosomeIndexStructure>();
		Path p = new Path(cramIndexFileName);
		Configuration conf = new Configuration();
		FileSystem fs = null;
		try {
			fs = p.getFileSystem(conf);
		} catch (IOException e1) {
			e1.printStackTrace();
		}

		try {
			if (!fs.exists(p)) {
				writeIndex(indexForChromosome(new Path(cramFileName)));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		try {
			FSDataInputStream reader = p.getFileSystem(conf).open(p);
			LineReader lineReader = new LineReader(reader, conf);
			Text line = new Text();
			while (lineReader.readLine(line) > 0) {
				if (line.getLength() == 0)
					continue;
				String[] str = line.toString().split("\t");
				chromosomeIndexMap.put(Integer.parseInt(str[0]),
						new ChromosomeIndexStructure(line.toString()));
			}

			lineReader.close();
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public long getStart(String chrName) {
		if (chromosomeIndexMap == null)
			loadChromosomeIndex();
		int seqId = samFileHeader.getSequenceIndex(chrName);
		if (chromosomeIndexMap.containsKey(seqId))
			return chromosomeIndexMap.get(seqId).getStart();
		else
			return -1;
	}

	public long getEnd(String chrName) {
		if (chromosomeIndexMap == null)
			loadChromosomeIndex();
		int seqId = samFileHeader.getSequenceIndex(chrName);
		if (chromosomeIndexMap.containsKey(seqId))
			return chromosomeIndexMap.get(seqId).getEnd();
		else
			return -1;
	}
}
