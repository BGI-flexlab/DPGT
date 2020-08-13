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
package org.bgi.flexlab.gaea.tools.recalibrator.covariate;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.tools.recalibrator.ReadCovariates;

import java.util.HashMap;
import java.util.List;

public class ReadGroupCovariate implements RequiredCovariate {
	private final HashMap<String, Integer> readGroupLookupTable = new HashMap<String, Integer>();
	private final HashMap<Integer, String> readGroupReverseLookupTable = new HashMap<Integer, String>();
	private int currentId = 0;

	@Override
	public void initialize(RecalibratorOptions option) {
	}

	public void initializeReadGroup(SAMFileHeader mFileHeader) {
		List<SAMReadGroupRecord> rg = mFileHeader.getReadGroups();
		for (SAMReadGroupRecord r : rg) {
			keyForReadGroup(readGroupValueFromRG(r));
		}
	}

	@Override
	public void recordValues(final GaeaSamRecord read, final ReadCovariates values) {
		final String readGroupId = readGroupValueFromRG(read.getReadGroup());
		final int key = keyForReadGroup(readGroupId);
		final int l = read.getReadLength();
		for (int i = 0; i < l; i++)
			values.addCovariate(key, key, key, i);
	}

	@Override
	public final Object getValue(final String str) {
		return str;
	}

	@Override
	public String formatKey(final int key) {
		return readGroupReverseLookupTable.get(key);
	}

	@Override
	public int keyFromValue(final Object value) {
		return keyForReadGroup((String) value);
	}

	private int keyForReadGroup(final String readGroupId) {
		if (!readGroupLookupTable.containsKey(readGroupId)) {
			readGroupLookupTable.put(readGroupId, currentId);
			readGroupReverseLookupTable.put(currentId, readGroupId);
			currentId++;
		}
		return readGroupLookupTable.get(readGroupId);
	}

	@Override
	public int maximumKeyValue() {
		return currentId - 1;
	}

	/**
	 * If the sample has a PU tag annotation, return that. If not, return the
	 * read group id.
	 */
	private String readGroupValueFromRG(final SAMReadGroupRecord rg) {
		final String platformUnit = rg.getPlatformUnit();
		return platformUnit == null ? rg.getId() : platformUnit;
	}
}
