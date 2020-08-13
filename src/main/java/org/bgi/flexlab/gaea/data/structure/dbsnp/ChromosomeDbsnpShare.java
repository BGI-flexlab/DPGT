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
package org.bgi.flexlab.gaea.data.structure.dbsnp;

import org.bgi.flexlab.gaea.data.structure.memoryshare.BioMemoryShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;

public class ChromosomeDbsnpShare extends BioMemoryShare {
	private final static int WINDOW_SIZE = VcfIndex.WINDOW_SIZE;
	private final static int CAPACITY = Long.SIZE / Byte.SIZE;

	public ChromosomeDbsnpShare() {
		super(1);
	}

	/*public long getStartPosition(int winNum, int winSize) {
		if (winNum * winSize >= getLength())
			throw new RuntimeException("position is more than chromosome length.");

		if (winSize == 0)
			winSize = WINDOW_SIZE;

		if (winSize % WINDOW_SIZE != 0)
			throw new RuntimeException("window size is not multiple for " + WINDOW_SIZE);

		int multipe = winSize / WINDOW_SIZE;

		int minWinNum = winNum * multipe;
		int maxWinNum = (winNum + 1) * multipe;
		
		if(minWinNum * CAPACITY >= fcSize)
			return -1;
		int end = maxWinNum * CAPACITY - 1;
		if(end >= fcSize)
			end = (fcSize - 1);
			
		byte[] indexs = getBytes(minWinNum * CAPACITY, end);

		long position = 0;

		for (int j = 0; j < indexs.length; j += CAPACITY) {
			for (int i = 0; i < CAPACITY; i++) {
				position <<= CAPACITY;
				position |= (indexs[j + i] & 0xff);
			}
			if (position != 0)
				return position;

			position = 0;
		}

		return -1;
	}*/
	
	public long getStartPosition(int winNum, int winSize) {
		return getStartPosition(winNum,winNum+1,winSize);
	}
	
	public long getStartPosition(int startWinNum, int endWinNum, int winSize) {
		if (startWinNum * winSize >= getLength())
			throw new RuntimeException("position is more than chromosome length.");

		if (winSize == 0)
			winSize = WINDOW_SIZE;

		if (winSize % WINDOW_SIZE != 0)
			throw new RuntimeException("window size is not multiple for " + WINDOW_SIZE);

		int multipe = winSize / WINDOW_SIZE;

		int minWinNum = startWinNum * multipe;
		int maxWinNum = endWinNum * multipe;
		
		if(minWinNum * CAPACITY >= fcSize)
			return -1;
		int end = maxWinNum * CAPACITY - 1;
		if(end >= fcSize)
			end = (fcSize - 1);
			
		byte[] indexs = getGA4GHBytes(minWinNum * CAPACITY, end);

		long position = 0;

		for (int j = 0; j < indexs.length; j += CAPACITY) {
			for (int i = 0; i < CAPACITY; i++) {
				position <<= CAPACITY;
				position |= (indexs[j + i] & 0xff);
			}
			if (position != 0)
				return position;

			position = 0;
		}

		return -1;
	}

	public long getStartPosition(int winNum) {
		return getStartPosition(winNum, WINDOW_SIZE);
	}
}
