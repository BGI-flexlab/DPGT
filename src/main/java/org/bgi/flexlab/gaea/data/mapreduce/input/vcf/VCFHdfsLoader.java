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
package org.bgi.flexlab.gaea.data.mapreduce.input.vcf;

import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.bgi.flexlab.gaea.data.mapreduce.util.HdfsFileManager;
import org.bgi.flexlab.gaea.data.structure.vcf.AbstractVCFLoader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import java.io.IOException;
import java.io.InputStream;

public class VCFHdfsLoader extends AbstractVCFLoader{

	private WrapSeekable seekableStream;
	
	public VCFHdfsLoader(String dbSNP) throws IllegalArgumentException, IOException {
		super(dbSNP);
	}
	
	@Override
	public void seek(long pos) {
		this.pos = pos;
		try {
			seekableStream.seek(pos);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}
		iterator = new AsciiLineReaderIterator(new AsciiLineReader(seekableStream));
	}

	@Override
	public void getSeekableStream() {
		// TODO Auto-generated method stub
		try {
			seekableStream = WrapSeekable.openPath(new Configuration(), new Path(super.input));
		} catch (IllegalArgumentException | IOException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}
	}

	@Override
	public InputStream getInputStream(String input) {
		// TODO Auto-generated method stub
		return HdfsFileManager.getInputStream(new Path(super.input), new Configuration());
	}
}
