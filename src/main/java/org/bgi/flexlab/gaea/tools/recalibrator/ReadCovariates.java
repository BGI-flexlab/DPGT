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
package org.bgi.flexlab.gaea.tools.recalibrator;

import org.bgi.flexlab.gaea.util.EventType;

public class ReadCovariates {
	private int[][][] mkeys = null;
	private boolean largeKeys = false;
	private int currCovIndex = 0;
	
	public ReadCovariates(int readLength, int numberOfCovariates,boolean largeKeys) {
		this.largeKeys = largeKeys;
		int length = 1;
		if(largeKeys)
			length = EventType.values().length;
		initialize(length,readLength,numberOfCovariates);
	}
	
	public ReadCovariates(int readLength, int numberOfCovariates) {
		this(readLength,numberOfCovariates,false);
	}
	
	private void initialize(int length,final int readLength, final int numberOfCovariates) {
		mkeys = new int[length][readLength][numberOfCovariates];
	}

	public void setCovariateIndex(final int index) {
		currCovIndex = index;
	}

	private void setCovariate(final int readOffset, final int... keys) {
		int len = keys.length;
		for (int i = 0; i < len; i++) {
			mkeys[i][readOffset][currCovIndex] = keys[i];
		}
	}

	public void addCovariate(final int mismatch, final int insertion, final int deletion, final int readOffset) {
		if (!largeKeys)
			setCovariate(readOffset, mismatch);
		else {
			setCovariate(readOffset, mismatch, insertion, deletion);
		}
	}

	public int[] getKeySet(final int readPosition, final EventType errorModel) {
		if (largeKeys)
			return mkeys[errorModel.index][readPosition];
		if (errorModel.index != EventType.SNP.index)
			throw new RuntimeException("model not match");
		return mkeys[0][readPosition];
	}

	public int[][] getKeySet(final EventType errorModel) {
		if (largeKeys)
			return mkeys[errorModel.index];
		if (errorModel.index != EventType.SNP.index)
			throw new RuntimeException("model not match");
		return mkeys[0];
	}

	public int[] getMismatchesKeySet(final int readPosition) {
		return mkeys[0][readPosition];
	}

	public int[] getInsertionsKeySet(final int readPosition) {
		if (!largeKeys)
			throw new RuntimeException("insertion model not match");
		return mkeys[EventType.Insertion.index][readPosition];
	}

	public int[] getDeletionsKeySet(final int readPosition) {
		if (!largeKeys)
			throw new RuntimeException("deletetion model not match");
		return mkeys[EventType.Deletion.index][readPosition];
	}

	public int[][] getMismatchesKeySet() {
		return mkeys[0];
	}

	public int[][] getInsertionsKeySet() {
		if (!largeKeys)
			throw new RuntimeException("insertion model not match");
		return mkeys[EventType.Insertion.index];
	}

	public int[][] getDeletionsKeySet() {
		if (!largeKeys)
			throw new RuntimeException("deletetion model not match");
		return mkeys[EventType.Deletion.index];
	}
	
	public int[][][] getMKeys(){
		return mkeys;
	}

	public void clear() {
		mkeys = null;
	}
}
