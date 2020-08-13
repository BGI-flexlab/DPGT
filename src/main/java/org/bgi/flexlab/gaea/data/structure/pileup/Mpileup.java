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
import org.bgi.flexlab.gaea.data.structure.pileup.filter.PileupFilter;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class Mpileup implements MpileupInterface<Pileup>{

	/**
	 * reads pool list
	 */
	private ReadsPool readsPool;

	/**
	 * sample -> pileup
	 */
	private Map<String, Pileup> pileups = new HashMap<>();

	/**
	 * end of pileup position
	 */
	private int end;

	/**
	 * position
	 */
	protected int position;

	/**
	 * tmp read
	 */
	private AlignmentsBasic tmpRead = null;

	/**
	 * make the pileup filter to static so that we do not need to store this in each pileup
	 */
	private PileupFilter filter;

	/**
	 * constructor
	 * @param readsPool reads
	 * @param position start
	 * @param end end
	 */
	public Mpileup(ReadsPool readsPool, int position, int end, PileupFilter filter) {
		this.readsPool = readsPool;
		this.position = position;
		this.end = end;
		this.filter = filter;
	}

	/**
	 * add reads
	 * @param read read
	 * @param pos position
	 */
	protected void addReads2Pileup(AlignmentsBasic read, int pos) {
		Pileup pileup = pileups.get(read.getSample());
		if (pileup == null) {
			pileup = new Pileup();
			pileup.setPosition(pos);
			pileups.put(read.getSample(), pileup);
		}
		pileup.addReads(read);
	}

	/**
	 *
	 * @return is Empty of all sample pileup
	 */
	public boolean allEmpty() {
		boolean allEmpty = true;
		for (String sample : pileups.keySet()) {
			Pileup pileup = pileups.get(sample);
			if (!pileup.isEmpty()) {
				allEmpty = false;
				break;
			}
		}

		return allEmpty;
	}

	/**
	 *
	 * @param minPosition
	 * @param size
	 * @return min position
	 */
	public int forwardPosition(int minPosition, int size) {
		int minimumPosition = Integer.MAX_VALUE;
		// forward position
		@SuppressWarnings("rawtypes")
		java.util.Iterator it = pileups.entrySet().iterator();
		while(it.hasNext()){
			@SuppressWarnings("rawtypes")
			java.util.Map.Entry entry = (java.util.Map.Entry) it.next();
			Pileup pileup = (Pileup)entry.getValue();
			if (pileup.getPosition() == minPosition) {
				pileup.forwardPosition(size);
			}

			if(pileup.getPosition() < minimumPosition && !pileup.isEmpty())
				minimumPosition = pileup.getPosition();
			else if(pileup.isEmpty()){
				it.remove();
			}
		}

		return minimumPosition;
	}

	/**
	 *
	 * @param minPosition
	 * @return
	 */
	public int addReads(int minPosition) {
		int currentPosition = minPosition;

		if (readsPool.hasNext() && tmpRead == null) {
			tmpRead = readsPool.next();
		}

		while (tmpRead != null) {
			AlignmentsBasic read = tmpRead;

			if (currentPosition == Integer.MAX_VALUE
					|| currentPosition == read.getPosition()) {
				if (currentPosition == Integer.MAX_VALUE)
					currentPosition = read.getPosition();
				addReads2Pileup(read, currentPosition);
				if (readsPool.hasNext())
					tmpRead = readsPool.next();
				else
					tmpRead = null;
			} else {
				break;
			}
		}
		return currentPosition;
	}

	/**
	 *
	 * @param minPosition
	 * @param posPlps
	 */
	public void syn(int minPosition,Map<String, Pileup> posPlps) {
		@SuppressWarnings("rawtypes")
		Iterator it = pileups.entrySet().iterator();
		while(it.hasNext()){
			@SuppressWarnings("rawtypes")
			java.util.Map.Entry entry = (java.util.Map.Entry) it.next();
			Pileup plp = (Pileup)entry.getValue();
			if (plp.getPosition() == minPosition) {
				plp.calculateBaseInfo();
				if(filter != null) {
					plp.setFilterPileup(filter.filter(plp));
				}
				String sample = (String)entry.getKey();
				posPlps.put(sample, plp);
			}
		}
	}


	public Map<String, Pileup> getNextPosPileup() {
		if (position > end)
			return null;
		//System.err.println("forward position:" + position);
		int minPosition = forwardPosition(position, 1);
		//System.err.println("add reads:" + minPosition);
		position = addReads(minPosition);

		if (minPosition != Integer.MAX_VALUE && position != minPosition)
			throw new RuntimeException("error in" + position + "\t" + minPosition);
		if (position > end || allEmpty())
			return null;
		//System.err.println("syn pileups:" + position);
		Map<String, Pileup> posPlps = new HashMap<String, Pileup>();
		syn(position, posPlps);

		return posPlps;
	}

	public int getMinPositionInPlp() {
		int minPosition = Integer.MAX_VALUE;

		for (String sample : pileups.keySet()) {
			Pileup pileup = pileups.get(sample);
			if (minPosition > pileup.getPosition() && !pileup.isEmpty()) {
				minPosition = pileup.getPosition();
			}
		}

		return minPosition;
	}

	public Pileup joinPileups() {
		Pileup joinedPileup = new Pileup();
		joinedPileup.setPosition(position);
		for(Pileup pileup : pileups.values()) {
			if(joinedPileup.getFilteredPileup() != null)
				joinedPileup.getFilteredPileup().addAll(pileup.getFilteredPileup());
			joinedPileup.getTotalPileup().addAll(pileup.getTotalPileup());
		}
		return joinedPileup;
	}

	public Map<String, Pileup> getCurrentPosPileup() {
		return pileups;
	}

	public int getPosition() {
		return position;
	}

	public int getSize() {
		return pileups.size();
	}
	
	public void clear(){
		pileups.clear();
	}

	public int totalDepth(boolean isFiltered) {
		int depth = 0;
		for(String sample : pileups.keySet()) {
			depth += pileups.get(sample).depthOfCoverage(isFiltered);
		}
		return depth;
	}
}

