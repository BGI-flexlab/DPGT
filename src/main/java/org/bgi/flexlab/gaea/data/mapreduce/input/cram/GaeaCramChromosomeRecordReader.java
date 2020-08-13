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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2010 Aalto University 
 *
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.mapreduce.input.cram;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;

import java.io.IOException;

public class GaeaCramChromosomeRecordReader extends GaeaCramRecordReader {
	private int sequenceId = Integer.MIN_VALUE;
	private int prevSeqId = -1;
	public final static String CHROMOSOME = "chromosome.name";

	public void initialize(InputSplit inputSplit, TaskAttemptContext context)
			throws IOException {
		super.initialize(inputSplit, context);

		FileSplit split = (FileSplit) inputSplit;
		final Path file = split.getPath();

		String chrName = context.getConfiguration().get(CHROMOSOME);
		String indexPath = context.getConfiguration().get("cram.index.path");

		if (chrName != null) {
			ChromosomeIndex chromosome = null;
			if (indexPath == null)
				chromosome = new ChromosomeIndex(file.toString());
			else
				chromosome = new ChromosomeIndex(file.toString(), indexPath
						+ "/" + file.getName() + ".crai");
			chromosome.setHeader(samFileHeader);
			start = chromosome.getStart(chrName);
			length = chromosome.getEnd(chrName) - start;
			
			sequenceId = samFileHeader.getSequenceIndex(chrName);
			seekableStream.seek(start);
		}
	}

	@Override
	public boolean nextKeyValue() {
		/* get new container */
		boolean res = super.nextKeyValue();
		int currSequenceId = record.get().getReferenceIndex();
		
		if(prevSeqId == -1){
			if(currSequenceId != sequenceId)
				throw new RuntimeException("Error index start!!");
			prevSeqId = currSequenceId;
		}else if(prevSeqId == sequenceId && currSequenceId != sequenceId){
			return false;
		}

		return res;
	}
}
