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
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.LineReader;

import java.io.*;

public class SamFileHeaderCodec {

	/**
	 * read sam header from file
	 * 
	 * @param LineReader for htsjdk
	 * @return SAMFileHeader
	 */
	public static SAMFileHeader readHeader(LineReader reader) {
		SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
		SAMFileHeader mHeader = codec.decode(reader, null);
		return mHeader;
	}
	
	/**
	 * write sam header to file for BufferedWriter
	 * 
	 * @param SAMFileHeader
	 * @param BufferedWriter
	 * @return
	 */
	public static void writeHeader(SAMFileHeader header, BufferedWriter writer)
			throws IOException {
		final Writer sw = new StringWriter();
		new SAMTextHeaderCodec().encode(sw, header);

		writer.write(sw.toString());
		writer.close();
	}

	/**
	 * write sam header to file for OutputStream
	 * 
	 * @param SAMFileHeader
	 * @param OutputStream
	 * @return
	 */
	public static void writeHeader(SAMFileHeader header, OutputStream out)
			throws IOException {
		writeHeader(header, new BufferedWriter(new OutputStreamWriter(out)));
	}
}
