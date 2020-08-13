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
package org.bgi.flexlab.gaea.data.structure.vcf;

import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.tribble.util.ParsingUtils;
import org.tukaani.xz.SeekableFileInputStream;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;

public class VCFLocalLoader extends AbstractVCFLoader{

	private SeekableFileInputStream seekableStream;

	private long pos;

	public VCFLocalLoader(String dbSNP) throws IOException {
		super(dbSNP);
	}

	
	@Override
	public void seek(long pos){
		this.pos = pos;
		try {
			seekableStream.seek(this.pos);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e.toString());
		}
		iterator = new AsciiLineReaderIterator(new AsciiLineReader(seekableStream));
	}


	@Override
	public void getSeekableStream() {
		// TODO Auto-generated method stub
		try {
			seekableStream = new SeekableFileInputStream(input);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e.toString());
		}
	}


	@Override
	public InputStream getInputStream(String input) throws IOException {
		// TODO Auto-generated method stub
		return ParsingUtils.openInputStream(input);
	}

}
