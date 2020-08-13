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
package org.bgi.flexlab.gaea.data.structure.pileup;

import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.util.MathUtils;

import java.util.ArrayList;

public class Pileup implements PileupInterface<PileupReadInfo> {

	public static int MAX_DEPTH = 50000;

	/**
	 * pileup struct;
	 */
	private ArrayList<PileupReadInfo> plp;

	/**
	 * filtered pileup
	 */
	private ArrayList<PileupReadInfo> filterPileup;

	/**
	 * position
	 */
	private int position;

	/**
	 * count of next base is deletion
	 */
	private int nextMatchCount;

	/**
	 * count of next base is insertion
	 */
	private int nextInsertionCount;

	/**
	 * count of next base is deletion
	 */
	private int nextDeletionCount;

	/**
	 * count of deletion base
	 */
	private int deletionCount;

	public Pileup() {
		position = -1;
		plp = new ArrayList<>();
		filterPileup = new ArrayList<>();
	}

	/**
	 * add readInfo to pileup
	 *
	 * @param readInfo
	 *            read info in AlignmentsBasic
	 */
	public void addReads(AlignmentsBasic readInfo) {
		PileupReadInfo read = new PileupReadInfo(readInfo);
		if (position >= read.getPosition() && position <= read.getEnd() && plp.size() < MAX_DEPTH) {
			plp.add(read);
		} else {
			if (plp.size() == 0) {
				position = read.getPosition();
				plp.add(read);
			} else if (position < read.getPosition() || position > read.getEnd()) {
				throw new RuntimeException("add read to plp error：" + read.getReadInfo().toString());
			}
		}
	}

	/**
	 * remove proccessed reads
	 */
	public void remove() {
		for (int i = 0; i < plp.size(); i++) {
			PileupReadInfo posRead = plp.get(i);

			if (position > posRead.getEnd()) {
				plp.remove(i);
				i--;
			}
		}
	}

	/**
	 * forward the pileup
	 * @param size forward loc length
	 */
	public void forwardPosition(int size) {
		position += size;

		remove();

		if (isEmpty())
			position = Integer.MAX_VALUE;
	}

	/**
	 * calculate read base info
	 */
	@Override
	public void calculateBaseInfo() {
		deletionCount = 0;
		nextDeletionCount = 0;
		nextInsertionCount = 0;
		if (position != Integer.MAX_VALUE) {
			filterPileup = new ArrayList<>();
			for (int i = 0; i < plp.size(); i++) {
				PileupReadInfo posRead = plp.get(i);
				posRead.calculateQueryPosition(position);
				if (posRead.isDeletionBase())
					deletionCount++;
				if (posRead.isNextDeletionBase())
					nextDeletionCount++;
				if (posRead.isNextInsertBase())
					nextInsertionCount++;
				if (posRead.isNextMatchBase())
					nextMatchCount++;
			}
		}
	}

	/**
	 * get pileup bases
	 * @return pileup bases
	 */
	public char[] getBases() {
		int i = 0;
		char[] bases = new char[plp.size()];
		for(PileupReadInfo read : plp) {
			bases[i++] = read.getBase();
		}

		return bases;
	}

	/**
	 * get pileup coverage
	 * FIXME:: maybe need to be fixed with filters.
	 * @return coverage depth
	 */
	public int depthOfCoverage(boolean isFiltered) {
		ArrayList<PileupReadInfo> pileupElements = plp;
		if(isFiltered)
			pileupElements = filterPileup;

		int depth = pileupElements.size();
		/*for(PileupReadInfo p : pileupElements) {
			depth += p.isDeletionBase() ?  : 1;
		}*/

		return depth;
	}

	public int getNumberOfDeletions() {
		return deletionCount;
	}

	/**
	 * is plp empty
	 *
	 * @return
	 */
	public boolean isEmpty() {
		return plp.size() == 0;
	}

	/**
	 * get position
	 *
	 * @return
	 */
	public int getPosition() {
		return position;
	}

	/**
	 * set position
	 *
	 * @param position
	 *            position
	 */
	public void setPosition(int position) {
		this.position = position;
	}

	/**
	 *
	 * @return plp array list.
	 */
	public ArrayList<PileupReadInfo> getTotalPileup() {
		return plp;
	}


	public ArrayList<PileupReadInfo> getFilteredPileup() {
		return filterPileup;
	}

	public void setFilterPileup(ArrayList<PileupReadInfo> filterPileup) {
		this.filterPileup = filterPileup;
	}

	public int getDeletionCount() {
		return deletionCount;
	}

	public double getDeletionRate() {
		return deletionCount / (double) depthOfCoverage(false);
	}

	public int getNextDeletionCount() {
		return nextDeletionCount;
	}

	public int getNextInsertionCount() {
		return nextInsertionCount;
	}

	public double getNextIndelRate() {
		return (nextDeletionCount + nextInsertionCount) / (double) nextMatchCount;
	}

	public int getNumberOfElements() {
		return plp.size();
	}
}
