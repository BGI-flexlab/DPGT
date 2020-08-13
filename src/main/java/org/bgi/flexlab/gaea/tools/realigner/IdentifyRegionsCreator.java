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
package org.bgi.flexlab.gaea.tools.realigner;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.pileup.Mpileup;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;
import org.bgi.flexlab.gaea.data.structure.pileup.ReadsPool;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RealignerOptions;
import org.bgi.flexlab.gaea.tools.realigner.event.Event;
import org.bgi.flexlab.gaea.tools.realigner.event.Event.EVENT_TYPE;
import org.bgi.flexlab.gaea.tools.realigner.event.EventPair;

import java.util.ArrayList;
import java.util.Map;

public class IdentifyRegionsCreator {
	private RealignerOptions option = null;
	private ArrayList<VariantContext> knowIndels = null;
	private ChromosomeInformationShare chr = null;
	private ArrayList<GaeaSamRecord> records = null;
	private GenomeLocationParser parser = null;
	private ArrayList<GenomeLocation> intervals = null;
	private int maxIntervalSize = 500;

	public IdentifyRegionsCreator(RealignerOptions option, ArrayList<GaeaSamRecord> records, SAMFileHeader mHeader,
			ChromosomeInformationShare chr, ArrayList<VariantContext> knowIndels) {
		this.records = records;
		this.knowIndels = knowIndels;
		this.parser = new GenomeLocationParser(mHeader.getSequenceDictionary());
		this.chr = chr;
		this.option = option;
		this.intervals = new ArrayList<GenomeLocation>();
		maxIntervalSize = option.getMaxInterval();
	}

	public Event getEvent(VariantState state, int chrIndex, Pileup pileup, int position) {
		boolean hasIndel = state.isIndel();
		boolean hasInsertion = state.isInsertion();
		boolean hasSNP = state.isSNP();
		int furthestPosition = state.getFurthestPosition();

		byte refBase = chr.getBinaryBase(position - 1);

		boolean lookForMismatchEntropy = option.getMismatchThreshold() > 0 && option.getMismatchThreshold() <= 1 ? true
				: false;

		long totalQuality = 0, mismatchQuality = 0;
		
		pileup.calculateBaseInfo();

		for (PileupReadInfo p : pileup.getTotalPileup()) {
			furthestPosition = Math.max(furthestPosition, p.getAlignmentEnd());

			if (p.isDeletionBase() || p.isNextInsertBase()) {
				hasIndel = true;
				if (p.isNextInsertBase())
					hasInsertion = true;
			} else if (lookForMismatchEntropy) {
				totalQuality += p.getBaseQuality();
				if (p.getBase() != refBase)
					mismatchQuality += p.getBaseQuality();
			}
		}

		if (lookForMismatchEntropy && pileup.getNumberOfElements() >= option.getMinReads()
				&& (double) mismatchQuality / (double) totalQuality >= option.getMismatchThreshold()){
			hasSNP = true;
		}

		if ((!hasIndel && !hasSNP) || (furthestPosition == -1)){
			return null;
		}

		GenomeLocation location = parser.createGenomeLocation(chrIndex, position);
		if (hasInsertion)
			location = parser.createGenomeLocation(location.getContig(), location.getStart(), location.getStart() + 1);

		EVENT_TYPE type = hasSNP ? (hasIndel ? EVENT_TYPE.BOTH : EVENT_TYPE.POINT_EVENT) : EVENT_TYPE.INDEL_EVENT;

		return new Event(location, furthestPosition, type);
	}

	public void setEventPair(EventPair pair, Event event) {
		if (event != null) {
			if (pair.getLeft() == null)
				pair.setLeft(event);
			else if (pair.getRight() == null) {
				if (pair.getLeft().canBeMerge(event)) {
					pair.getLeft().merge(event, option.getSNPWindowSize());
				} else {
					pair.setRight(event);
				}
			} else {
				if (pair.getRight().canBeMerge(event)) {
					pair.getRight().merge(event, option.getSNPWindowSize());
				} else {
					if (pair.getRight().isValidEvent(parser, maxIntervalSize)) {
						pair.getIntervals().add(pair.getRight().getLocation(parser));
					}
					pair.setRight(event);
				}
			}
		}
	}

	public void setIntervals(EventPair pair) {
		if (pair.getLeft() != null && pair.getLeft().isValidEvent(parser, maxIntervalSize))
			pair.getIntervals().add(pair.getLeft().getLocation(parser));
		if (pair.getRight() != null && pair.getRight().isValidEvent(parser, maxIntervalSize))
			pair.getIntervals().add(pair.getRight().getLocation(parser));

		for (GenomeLocation location : pair.getIntervals())
			intervals.add(location);
	}

	public ArrayList<GenomeLocation> getIntervals() {
		return intervals;
	}

	public void regionCreator(int chrIndex,int start, int end) {
		if(records .size() == 0)
			return;
		
		EventPair pair = new EventPair(null, null);
		ReadsPool pool = new ReadsPool(records.iterator(), null);
		
		int regionStart = records.get(0).getAlignmentStart();
		Mpileup mpileup = new Mpileup(pool, regionStart, end-1,null);

		Map<String, Pileup> pileups = mpileup.getNextPosPileup();
		
		if(pileups == null)
			return;

		if (pileups.size() != 1) {
			throw new RuntimeException("realigner pileup is more than one sample?");
		}

		while (pileups != null) {
			int currPosition = mpileup.getPosition()+1;
			for (String key : pileups.keySet()) {
				VariantState state = new VariantState();
				state.filterVariant(knowIndels, currPosition);
				Event event = getEvent(state,chrIndex, pileups.get(key), currPosition);
				setEventPair(pair, event);
			}
			pileups = mpileup.getNextPosPileup();
		}

		mpileup.clear();

		setIntervals(pair);
	}
}

