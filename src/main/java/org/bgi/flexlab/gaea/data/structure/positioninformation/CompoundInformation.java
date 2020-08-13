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
package org.bgi.flexlab.gaea.data.structure.positioninformation;

import org.bgi.flexlab.gaea.data.structure.bam.SAMInformationBasic;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.util.CigarState;
import org.bgi.flexlab.gaea.util.SamRecordDatum;
import org.bgi.flexlab.gaea.util.SystemConfiguration;

/**
 * a base class for compound information, which includes window information 
 * and position information (so far)
 * @param <T extends SAMInformationBasic>
 */
public class CompoundInformation<T extends SAMInformationBasic> {
	private T readInfo;
	
    private int refPosition = 0;
	
	private int readPosition = 0;
	
	private ChromosomeInformationShare chrInfo = null;
	
	private int windowStart = 0;
	
	private int windowSize = 0;
	
	/**
	 * Constructor for position information
	 */
	public CompoundInformation( T readInfo, ChromosomeInformationShare chrInfo,
			int windowStart, int refPosition, int readPosition) {
		initialize(windowStart, readInfo, chrInfo);
		this.refPosition = refPosition;
		this.readPosition = readPosition;
	}
	
	/**
	 * Constructor for window information
	 */
	public CompoundInformation(int windowStart, int winSize, T readInfo, ChromosomeInformationShare chrInfo) {
		initialize(windowStart, readInfo, chrInfo);
		this.windowSize = winSize;
	}
	
	private void initialize(int windowStart, T readInfo, ChromosomeInformationShare chrInfo) {
		this.windowStart = windowStart;
		this.readInfo = readInfo;
		this.chrInfo = chrInfo;
	}
	
	public boolean eligiblePos() {
		return (readInfo.getFlag() & (SystemConfiguration.BAM_FQCFAIL | 
				SystemConfiguration.BAM_FSECONDARY | SystemConfiguration.BAM_FDUP)) == 0 ;
	}

	public T getReadInfo() {
		return readInfo;
	}

	public int getRgIndex() {
		return ((SamRecordDatum) readInfo).getRgIndex();
	}
	
	public int getMappingQual() {
		return readInfo.getMappingQual();
	}
	
	public CigarState getCigarState() {
		return readInfo.getCigarState();	
	}
	
	public boolean isDup() {
		return readInfo.isDup();
	}

	public ChromosomeInformationShare getChrInfo() {
		return chrInfo;
	}

	public int getWindowStart() {
		return windowStart;
	}
	
	public int distBetweenRefPosAndWinStart() {
		return refPosition - windowStart;
	}
	
	public byte getBinaryBase() {
		return readInfo.getBinaryBase(readPosition);
	}
	
	public int getRefPosition() {
		return refPosition;
	}

	public int getReadPosition() {
		return readPosition;
	}

	public int getWindowSize(){
		return windowSize;
	}
}
