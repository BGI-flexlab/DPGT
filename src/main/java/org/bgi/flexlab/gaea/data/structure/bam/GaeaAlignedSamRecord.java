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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class GaeaAlignedSamRecord {
	public static final String ORIGINAL_CIGAR_TAG = "OC";
	public static final String ORIGINAL_POSITION_TAG = "OP";
	public int MAX_POSITION_MOVE_ALLOWED = 200;
	private GaeaSamRecord read;
	private byte[] readBases = null;
	private byte[] baseQualities = null;
	private Cigar newCigar = null;
	private int newStart = -1;
	private int mismatchScore = 0;
	private int alignerMismatchScore = 0;
	private boolean ignoredOriginTag = false;

	public GaeaAlignedSamRecord(GaeaSamRecord read) {
		this.read = read;
		this.mismatchScore = 0;
	}

	public GaeaAlignedSamRecord(GaeaSamRecord read, int max_allow) {
		this(read);
		this.MAX_POSITION_MOVE_ALLOWED = max_allow;
	}

	public GaeaAlignedSamRecord(GaeaSamRecord read, int max_allow, boolean noTag) {
		this(read, max_allow);
		this.ignoredOriginTag = noTag;
	}

	public void setMismatchScore(int score) {
		this.mismatchScore = score;
	}

	public void setNewStart(int start) {
		this.newStart = start;
	}

	public void setAlignerMismatchScore(int score) {
		this.alignerMismatchScore = score;
	}

	public void setCigar(Cigar cigar) {
		setCigar(cigar, true);
	}

	private void setCigar(Cigar cigar, boolean reclip) {
		if (cigar == null) {
			return;
		}

		if(reclip && hasClipOperator(read.getCigar())){
			cigar = GaeaCigar.reclipCigar(cigar, read);
		}

		if (cigar.equals(read.getCigar())) {
			cigar = null;
		}

		newCigar = cigar;
	}
	
	private boolean hasClipOperator(Cigar cigar){
		for(CigarElement ce : cigar.getCigarElements()){
			if(ce.getOperator() == CigarOperator.S)
				return true;
		}
		
		return false;
	}

	public boolean statusFinalize() {
		if (newCigar == null) {
			return false;
		}

		newStart = getAlignmentStart();
		if (Math.abs(newStart - read.getAlignmentStart()) > this.MAX_POSITION_MOVE_ALLOWED){
			return false;
		}

		if (!ignoredOriginTag) {
			read.setAttribute(ORIGINAL_CIGAR_TAG, read.getCigar().toString());
			if (newStart != read.getAlignmentStart())
				read.setAttribute(ORIGINAL_POSITION_TAG,
						read.getAlignmentStart());
		}

		read.setCigar(newCigar);
		read.setAlignmentStart(newStart);

		return true;
	}

	public byte[] getReadBases() {
		if (readBases == null)
			getUnclippedInformation();
		return readBases;
	}

	public byte[] getReadQualities() {
		if (baseQualities == null)
			getUnclippedInformation();
		return baseQualities;
	}

	private void getUnclippedInformation() {
		readBases = new byte[read.getReadLength()];
		baseQualities = new byte[read.getReadLength()];

		int startIndex = 0, baseCount = 0;

		byte[] reads = read.getReadBases();
		byte[] qualities = read.getBaseQualities();

		for (CigarElement element : read.getCigar().getCigarElements()) {
			int eLength = element.getLength();
			switch (element.getOperator()) {
			case S:
				startIndex += eLength;
				break;
			case M:
			case X:
			case EQ:
			case I:
				System.arraycopy(reads, startIndex, readBases, baseCount,
						eLength);
				System.arraycopy(qualities, startIndex, baseQualities,
						baseCount, eLength);
				startIndex += eLength;
				baseCount += eLength;
			default:
				break;
			}
		}

		if (startIndex != baseCount) {
			byte[] baseTemp = new byte[baseCount];
			byte[] qualTemp = new byte[baseCount];
			System.arraycopy(readBases, 0, baseTemp, 0, baseCount);
			System.arraycopy(baseQualities, 0, qualTemp, 0, baseCount);
			readBases = baseTemp;
			baseQualities = qualTemp;
		}
	}

	public GaeaSamRecord getRead() {
		return this.read;
	}

	public Cigar getCigar() {
		return newCigar == null ? read.getCigar() : newCigar;
	}
	
	public Cigar getNewCigar(){
		return newCigar;
	}

	public int getAlignmentStart() {
		return newStart == -1 ? read.getAlignmentStart() : newStart;
	}

	public int getMismatchScore() {
		return this.mismatchScore;
	}

	public int getAlignerMismatchScore() {
		return this.alignerMismatchScore;
	}

	public int getReadLength() {
		return getReadBases().length;
	}
}
