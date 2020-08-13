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

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.sequenceplatform.NGSPlatform;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.tools.recalibrator.ReadCovariates;
import org.bgi.flexlab.gaea.util.BaseUtils;

import java.util.EnumSet;

public class CycleCovariate implements OptionalCovariate {
	private static final int MAXIMUM_CYCLE_VALUE = 1000;
	private static final int CUSHION_FOR_INDELS = 4;
	private String defaultPlatform = null;

	private static final EnumSet<NGSPlatform> DISCRETE_CYCLE_PLATFORMS = EnumSet.of(NGSPlatform.ILLUMINA,
			NGSPlatform.SOLID, NGSPlatform.PACBIO, NGSPlatform.COMPLETE_GENOMICS);
	private static final EnumSet<NGSPlatform> FLOW_CYCLE_PLATFORMS = EnumSet.of(NGSPlatform.LS454,
			NGSPlatform.ION_TORRENT);

	@Override
	public void initialize(RecalibratorOptions option) {
		if (option.DEFAULT_PLATFORM != null && !NGSPlatform.isKnown(option.DEFAULT_PLATFORM))
			throw new UserException("unknow platform for " + option.DEFAULT_PLATFORM);

		if (option.DEFAULT_PLATFORM != null)
			defaultPlatform = option.DEFAULT_PLATFORM;
	}

	@Override
	public void recordValues(GaeaSamRecord read, ReadCovariates values) {
		final int readLength = read.getReadLength();
		final NGSPlatform ngsPlatform = defaultPlatform == null ? read.getNGSPlatform()
				: NGSPlatform.fromReadGroupPL(defaultPlatform);

		if (DISCRETE_CYCLE_PLATFORMS.contains(ngsPlatform)) {
			final int readOrderFactor = read.getReadPairedFlag() && read.getSecondOfPairFlag() ? -1 : 1;
			final int increment;
			int cycle;
			if (read.getReadNegativeStrandFlag()) {
				cycle = readLength * readOrderFactor;
				increment = -1 * readOrderFactor;
			} else {
				cycle = readOrderFactor;
				increment = readOrderFactor;
			}

			final int MAX_CYCLE_FOR_INDELS = readLength - CUSHION_FOR_INDELS - 1;
			for (int i = 0; i < readLength; i++) {
				final int substitutionKey = keyFromCycle(cycle);
				final int indelKey = (i < CUSHION_FOR_INDELS || i > MAX_CYCLE_FOR_INDELS) ? -1 : substitutionKey;
				values.addCovariate(substitutionKey, indelKey, indelKey, i);
				cycle += increment;
			}
		} else if (FLOW_CYCLE_PLATFORMS.contains(ngsPlatform)) {

			final byte[] bases = read.getReadBases();

			final boolean multiplyByNegative1 = read.getReadPairedFlag() && read.getSecondOfPairFlag();

			int cycle = multiplyByNegative1 ? -1 : 1;

			if (!read.getReadNegativeStrandFlag()) { // Forward direction
				int iii = 0;
				while (iii < readLength) {
					while (iii < readLength && bases[iii] == (byte) 'T') {
						final int key = keyFromCycle(cycle);
						values.addCovariate(key, key, key, iii);
						iii++;
					}
					while (iii < readLength && bases[iii] == (byte) 'A') {
						final int key = keyFromCycle(cycle);
						values.addCovariate(key, key, key, iii);
						iii++;
					}
					while (iii < readLength && bases[iii] == (byte) 'C') {
						final int key = keyFromCycle(cycle);
						values.addCovariate(key, key, key, iii);
						iii++;
					}
					while (iii < readLength && bases[iii] == (byte) 'G') {
						final int key = keyFromCycle(cycle);
						values.addCovariate(key, key, key, iii);
						iii++;
					}
					if (iii < readLength) {
						if (multiplyByNegative1)
							cycle--;
						else
							cycle++;
					}
					if (iii < readLength && !BaseUtils.isRegularBase(bases[iii])) {
						final int key = keyFromCycle(cycle);
						values.addCovariate(key, key, key, iii);
						iii++;
					}

				}
			} else { // Negative direction
				int iii = readLength - 1;
				while (iii >= 0) {
					while (iii >= 0 && bases[iii] == (byte) 'T') {
						final int key = keyFromCycle(cycle);
						values.addCovariate(key, key, key, iii);
						iii--;
					}
					while (iii >= 0 && bases[iii] == (byte) 'A') {
						final int key = keyFromCycle(cycle);
						values.addCovariate(key, key, key, iii);
						iii--;
					}
					while (iii >= 0 && bases[iii] == (byte) 'C') {
						final int key = keyFromCycle(cycle);
						values.addCovariate(key, key, key, iii);
						iii--;
					}
					while (iii >= 0 && bases[iii] == (byte) 'G') {
						final int key = keyFromCycle(cycle);
						values.addCovariate(key, key, key, iii);
						iii--;
					}
					if (iii >= 0) {
						if (multiplyByNegative1)
							cycle--;
						else
							cycle++;
					}
					if (iii >= 0 && !BaseUtils.isRegularBase(bases[iii])) {
						final int key = keyFromCycle(cycle);
						values.addCovariate(key, key, key, iii);
						iii--;
					}
				}
			}
		} else {
			throw new UserException("The platform (" + read.getReadGroup().getPlatform()
					+ ") associated with read group " + read.getReadGroup()
					+ " is not a recognized platform. Implemented options are e.g. illumina, 454, and solid");
		}
	}

	@Override
	public Object getValue(String str) {
		return Integer.parseInt(str);
	}

	@Override
	public String formatKey(int key) {
		int cycle = key >> 1; // shift so we can remove the "sign" bit
		if ((key & 1) != 0) // is the last bit set?
			cycle *= -1; // then the cycle is negative
		return String.format("%d", cycle);
	}

	@Override
	public int keyFromValue(Object value) {
		return (value instanceof String) ? keyFromCycle(Integer.parseInt((String) value))
				: keyFromCycle((Integer) value);
	}

	@Override
	public int maximumKeyValue() {
		return (MAXIMUM_CYCLE_VALUE << 1) + 1;
	}

	private int keyFromCycle(final int cycle) {
		int result = Math.abs(cycle);
		result = result << 1;
		if (cycle < 0)
			result++;
		return result;
	}
}
