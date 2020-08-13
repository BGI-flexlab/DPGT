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
 * Copyright (c) 2009-2012 The Broad Institute
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
package org.bgi.flexlab.gaea.data.structure.bam;

import htsjdk.samtools.util.StringUtil;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;

import java.util.ArrayList;

public class GaeaSamRecordBin {
	private ArrayList<GaeaSamRecord> records = null;
	private byte[] refBases = null;
	private GenomeLocation location = null;
	private int REFERENCE_EXTEND = 30;
	private GenomeLocationParser parser = null;
	
	public GaeaSamRecordBin(GenomeLocationParser parser,int extend) {
		setParser(parser);
		setExtend(extend);
	}

	public GaeaSamRecordBin(GenomeLocationParser parser) {
		setParser(parser);
	}

	public void setParser(GenomeLocationParser parser) {
		this.parser = parser;
	}
	
	public void setExtend(int extend){
		this.REFERENCE_EXTEND = extend;
	}

	public void add(GaeaSamRecord read) {
		GenomeLocation loc = parser.createGenomeLocation(read);

		if (location == null)
			location = loc;
		else if (loc.getStop() > location.getStop()) {
			location = parser.createGenomeLocation(location.getContig(),
					location.getStart(), loc.getStop());
		}

		if(records == null)
			records = new ArrayList<GaeaSamRecord>();
		records.add(read);
	}

	public GenomeLocation getLocation() {
		return this.location;
	}

	public byte[] getReference(ChromosomeInformationShare chr,
			GenomeLocationParser parser) {
		if (refBases == null) {
			int start = location.getStart() > REFERENCE_EXTEND ? (location.getStart() - REFERENCE_EXTEND) : 1;
			int end = location.getStop() + REFERENCE_EXTEND > chr.getLength() ? chr.getLength() : (location.getStop() + REFERENCE_EXTEND);
			location = parser.createGenomeLocation(location.getContig(), start,
					end);
			refBases = chr.getGA4GHBaseBytes(start - 1, end - 1);

			StringUtil.toUpperCase(refBases);
		}
		return this.refBases;
	}

	public ArrayList<GaeaSamRecord> getReads() {
		return this.records;
	}

	public int size() {
		if(records == null)
			return 0;
		return records.size();
	}

	public void clear() {
		if(records != null)
			records.clear();
		refBases = null;
		location = null;
	}
}
