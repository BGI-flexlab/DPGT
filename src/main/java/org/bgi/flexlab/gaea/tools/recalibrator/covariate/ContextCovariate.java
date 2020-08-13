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
import org.bgi.flexlab.gaea.data.structure.bam.clipper.ReadClipper;
import org.bgi.flexlab.gaea.data.structure.bam.clipper.algorithm.ClippingWriteNs;
import org.bgi.flexlab.gaea.data.structure.bam.clipper.algorithm.ReadClippingAlgorithm;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.tools.recalibrator.ReadCovariates;
import org.bgi.flexlab.gaea.util.BaseUtils;

import java.util.ArrayList;

public class ContextCovariate implements OptionalCovariate {
	private static final int LENGTH_BITS = 4;
	private static final int LENGTH_MASK = (1 << LENGTH_BITS) - 1;
	private static final int MAX_DNA_CONTEXT = 13;

	private int mContextSize;
	private int iContextSize;
	private int mismatchesMask;
	private int indelsMask;

	private byte lowQuality;
	private ReadClippingAlgorithm algorithm;

	@Override
	public void initialize(RecalibratorOptions option) {
		mContextSize = option.MISMATCHES_CONTEXT_SIZE;
		iContextSize = option.INDELS_CONTEXT_SIZE;

		if (mContextSize <= 0 || mContextSize > MAX_DNA_CONTEXT)
			throw new UserException.BadArgumentValueException(mContextSize, 0, MAX_DNA_CONTEXT);

		if (iContextSize <= 0 || iContextSize > MAX_DNA_CONTEXT)
			throw new UserException.BadArgumentValueException(iContextSize, 0, MAX_DNA_CONTEXT);

		lowQuality = option.LOW_QUALITY_TAIL;

		mismatchesMask = createMask(mContextSize);
		indelsMask = createMask(iContextSize);

		algorithm = new ClippingWriteNs();
	}

	private int createMask(final int size) {
		int mask = 0;

		for (int i = 0; i < size; i++)
			mask = (mask << 2) | 3;

		return mask << LENGTH_BITS;
	}

	private ArrayList<Integer> contextWith(byte[] bases, int size, int mask) {
		final int readLength = bases.length;
		final ArrayList<Integer> keys = new ArrayList<Integer>(readLength);

		for (int i = 1; i < size && i <= readLength; i++)
			keys.add(-1);

		if (readLength < size)
			return keys;

		final int newBaseOffset = 2 * (size - 1) + LENGTH_BITS;

		int currentKey = keyFromContext(bases, 0, size);
		keys.add(currentKey);

		int currentNPenalty = 0;
		if (currentKey == -1) {
			currentKey = 0;
			currentNPenalty = size - 1;
			int offset = newBaseOffset;
			while (bases[currentNPenalty] != 'N') {
				final int baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentNPenalty]);
				currentKey |= (baseIndex << offset);
				offset -= 2;
				currentNPenalty--;
			}
		}

		for (int currentIndex = size; currentIndex < readLength; currentIndex++) {
			final int baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentIndex]);
			if (baseIndex == -1) {
				currentNPenalty = size;
				currentKey = 0;
			} else {
				currentKey = (currentKey >> 2) & mask;
				currentKey |= (baseIndex << newBaseOffset);
				currentKey |= size;
			}

			if (currentNPenalty == 0) {
				keys.add(currentKey);
			} else {
				currentNPenalty--;
				keys.add(-1);
			}
		}

		return keys;
	}

	@Override
	public void recordValues(GaeaSamRecord read, ReadCovariates values) {
		final byte[] originalBases = read.getReadBases().clone();
		ReadClipper clipper = new ReadClipper(read);
		final GaeaSamRecord clippedRead = clipper.clipLowQualityEnds(lowQuality, algorithm);

		final boolean negativeStrand = clippedRead.getReadNegativeStrandFlag();
		byte[] bases = clippedRead.getReadBases();
		if (negativeStrand)
			bases = BaseUtils.simpleReverseComplement(bases);

		ArrayList<Integer> mismatchKeys = contextWith(bases, mContextSize, mismatchesMask);
		ArrayList<Integer> indelKeys = contextWith(bases, iContextSize, indelsMask);

		final int readLength = bases.length;
		for (int i = 0; i < readLength; i++) {
			int indelKey = indelKeys.get(i);
			values.addCovariate(mismatchKeys.get(i), indelKey, indelKey, (negativeStrand ? readLength - i - 1 : i));
		}

		mismatchKeys.clear();
		indelKeys.clear();

		read.setReadBases(originalBases);
	}

	@Override
	public Object getValue(String str) {
		return str;
	}

	private static int keyFromContext(final byte[] dna, final int start, final int end) {
		int key = end - start;
		int bitOffset = LENGTH_BITS;
		for (int i = start; i < end; i++) {
			final int baseIndex = BaseUtils.simpleBaseToBaseIndex(dna[i]);
			if (baseIndex == -1)
				return -1;
			key |= (baseIndex << bitOffset);
			bitOffset += 2;
		}
		return key;
	}

	public static String contextFromKey(final int key) {
		if (key < 0)
			throw new UserException("dna conversion cannot handle negative numbers. Possible overflow?");

		int length = key & LENGTH_MASK;
		int mask = 3 << LENGTH_BITS;
		int offset = LENGTH_BITS;

		StringBuilder dna = new StringBuilder();
		for (int i = 0; i < length; i++) {
			final int baseIndex = (key & mask) >> offset;
			dna.append((char) BaseUtils.baseIndexToSimpleBase(baseIndex));
			mask = mask << 2;
			offset += 2;
		}

		return dna.toString();
	}

	@Override
	public String formatKey(int key) {
		if (key == -1)
			return null;

		return contextFromKey(key);
	}

	@Override
	public int keyFromValue(Object value) {
		return keyFromContext((String) value);
	}

	private int keyFromContext(String str) {
		return keyFromContext(str.getBytes(), 0, str.length());
	}

	@Override
	public int maximumKeyValue() {
		int length = Math.max(mContextSize, iContextSize);
		int key = 0;
		for (int i = 0; i < length; i++) {
			key <<= 2;
			key |= 3;
		}
		key <<= LENGTH_BITS;
		key |= length;

		return key;
	}
}
