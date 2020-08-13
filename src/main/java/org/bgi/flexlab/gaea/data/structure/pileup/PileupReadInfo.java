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
package org.bgi.flexlab.gaea.data.structure.pileup;


import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.bam.SAMCompressionInformationBasic;
import org.bgi.flexlab.gaea.util.BaseUtils;
import org.bgi.flexlab.gaea.util.CigarState;
import org.bgi.flexlab.gaea.util.SystemConfiguration;

public class PileupReadInfo {
	public static final byte DELETION_BASE = BaseUtils.D;
	public static final byte DELETION_QUAL = (byte) 16;
	public static final byte A_FOLLOWED_BY_INSERTION_BASE = (byte) 87;
	public static final byte C_FOLLOWED_BY_INSERTION_BASE = (byte) 88;
	public static final byte T_FOLLOWED_BY_INSERTION_BASE = (byte) 89;
	public static final byte G_FOLLOWED_BY_INSERTION_BASE = (byte) 90;
	/**
	 * read information
	 */
	protected AlignmentsBasic readInfo;

	/**
	 * position on reads
	 */
	protected int qpos;

	/**
	 * cigar state for cigar analysis
	 */
	protected CigarState cigarState = new CigarState();

	/**
	 * end for reads
	 */
	protected int end;

	/**
	 * sample info
	 */
	protected String sample;


	/**
	 * // what is the length of the event (insertion or deletion) *after* this base
	 */
	protected int eventLength = -1;


	/**
	 * if it is a deletion, we do not have information about the actual deleted bases in the read itself, so we fill the string with D's; for insertions we keep actual inserted bases
	 */
	protected String eventBases = null;


	/**
	 * constructor
	 */
	public PileupReadInfo() {
		readInfo = null;
		qpos = Integer.MIN_VALUE;
		end = readInfo.getPosition();
		sample = null;
	}

	/**
	 * constructor
	 * @param readInfo
	 */
	public PileupReadInfo(AlignmentsBasic readInfo) {
		this.readInfo = readInfo;
		cigarState.parseCigar(readInfo.getCigars());
		end = readInfo.calculateReadEnd();
		sample = readInfo.getSample();
	}

	/**
	 * calculate query position
	 * @param refPos
	 */
	public void calculateQueryPosition(int refPos) {
		qpos = cigarState.resolveCigar(refPos, readInfo.getPosition());

		eventLength = cigarState.getEventLength();
		if (isNextInsertBase())
			eventBases = readInfo.getReadBasesString(qpos + 1, eventLength);
		else
			eventBases = null;                  // ignore argument in any other case
	}

	public int calcAlignmentByteArrayOffset(final int alignmentStart, final int refLocus) {
		return calcAlignmentByteArrayOffset(isInsertionAtBeginningOfRead(), cigarState.isDeletionBase(), alignmentStart, refLocus);
	}

	public int calcAlignmentByteArrayOffset(final boolean isInsertionAtBeginningOfRead, final boolean isDeletion,
											final int alignmentStart, final int refLocus) {
		int pileupOffset = qpos;
		int[] cigars = readInfo.getCigars();

		// Special case for reads starting with insertion
		if (isInsertionAtBeginningOfRead)
			return 0;

		// Reassign the offset if we are in the middle of a deletion because of the modified representation of the read bases
		if (isDeletion) {
			pileupOffset = refLocus - alignmentStart;
			final int cigar = cigars[0];
			int cigarOp = (cigar & SystemConfiguration.BAM_CIGAR_MASK);
			int cigarLength = (cigar >> SystemConfiguration.BAM_CIGAR_SHIFT);
			if (cigarOp == SystemConfiguration.BAM_CSOFT_CLIP) {
				pileupOffset += cigarLength;
			}
		}

		int pos = 0;
		int alignmentPos = 0;

		for (int iii = 0; iii < cigars.length; iii++) {
			int cigar = cigars[iii];
			int cigarOp = (cigar & SystemConfiguration.BAM_CIGAR_MASK);
			int cigarLength = (cigar >> SystemConfiguration.BAM_CIGAR_SHIFT);

			switch (cigarOp) {
				case SystemConfiguration.BAM_CINS:
				case SystemConfiguration.BAM_CSOFT_CLIP:
					pos += cigarLength;
					if (pos >= pileupOffset) {
						return alignmentPos;
					}
					break;
				case SystemConfiguration.BAM_CDEL:
					if (!isDeletion) {
						alignmentPos += cigarLength;
					} else {
						if (pos + cigarLength - 1 >= pileupOffset) {
							return alignmentPos + (pileupOffset - pos);
						} else {
							pos += cigarLength;
							alignmentPos += cigarLength;
						}
					}
					break;
				case SystemConfiguration.BAM_CMATCH:
				case SystemConfiguration.BAM_CEQUAL:
				case SystemConfiguration.BAM_CDIFF:
					if (pos + cigarLength - 1 >= pileupOffset) {
						return alignmentPos + (pileupOffset - pos);
					} else {
						pos += cigarLength;
						alignmentPos += cigarLength;
					}
					break;
				case SystemConfiguration.BAM_CHARD_CLIP:
				case SystemConfiguration.BAM_CPAD:
				case SystemConfiguration.BAM_CREF_SKIP:
					break;
				default:
					throw new UserException("Unsupported cigar operator: " + cigarOp);
			}
		}

		return alignmentPos;
	}

	public static byte[] readToAlignmentByteArray(int[] cigars, final byte[] read) {

		int alignmentLength = 0;
		for (int iii = 0; iii < cigars.length; iii++) {

			int cigar = cigars[iii];
			int cigarOp = (cigar & SystemConfiguration.BAM_CIGAR_MASK);
			int cigarLength = (cigar >> SystemConfiguration.BAM_CIGAR_SHIFT);

			switch (cigarOp) {
				case SystemConfiguration.BAM_CDEL:
				case SystemConfiguration.BAM_CREF_SKIP:
				case SystemConfiguration.BAM_CMATCH:
				case SystemConfiguration.BAM_CEQUAL:
				case SystemConfiguration.BAM_CDIFF:
					alignmentLength += cigarLength;
					break;
				case SystemConfiguration.BAM_CINS:
				case SystemConfiguration.BAM_CSOFT_CLIP:
				case SystemConfiguration.BAM_CHARD_CLIP:
				case SystemConfiguration.BAM_CPAD:
					break;
				default:
					throw new UserException("Unsupported cigar operator: " + cigarOp);
			}
		}

		final byte[] alignment = new byte[alignmentLength];
		int alignPos = 0;
		int readPos = 0;
		for (int iii = 0; iii < cigars.length; iii++) {

			int cigar = cigars[iii];
			int cigarOp = (cigar & SystemConfiguration.BAM_CIGAR_MASK);
			int cigarLength = (cigar >> SystemConfiguration.BAM_CIGAR_SHIFT);

			switch (cigarOp) {
				case SystemConfiguration.BAM_CINS:
					if (alignPos > 0) {
						if (alignment[alignPos - 1] == BaseUtils.A) {
							alignment[alignPos - 1] = PileupReadInfo.A_FOLLOWED_BY_INSERTION_BASE;
						} else if (alignment[alignPos - 1] == BaseUtils.C) {
							alignment[alignPos - 1] = PileupReadInfo.C_FOLLOWED_BY_INSERTION_BASE;
						} else if (alignment[alignPos - 1] == BaseUtils.T) {
							alignment[alignPos - 1] = PileupReadInfo.T_FOLLOWED_BY_INSERTION_BASE;
						} else if (alignment[alignPos - 1] == BaseUtils.G) {
							alignment[alignPos - 1] = PileupReadInfo.G_FOLLOWED_BY_INSERTION_BASE;
						}
					}
				case SystemConfiguration.BAM_CSOFT_CLIP:
					for (int jjj = 0; jjj < cigarLength; jjj++) {
						readPos++;
					}
					break;
				case SystemConfiguration.BAM_CDEL:
				case SystemConfiguration.BAM_CREF_SKIP:
					for (int jjj = 0; jjj < cigarLength; jjj++) {
						alignment[alignPos] = PileupReadInfo.DELETION_BASE;
						alignPos++;
					}
					break;
				case SystemConfiguration.BAM_CMATCH:
				case SystemConfiguration.BAM_CEQUAL:
				case SystemConfiguration.BAM_CDIFF:
					for (int jjj = 0; jjj < cigarLength; jjj++) {
						try {
                            alignment[alignPos] = read[readPos];
                        } catch(ArrayIndexOutOfBoundsException e) {
						    throw  new RuntimeException(SAMCompressionInformationBasic.CigarToString(cigars) + "\t" +
                                    e.getMessage() + "\t" + alignPos + "\t" + readPos + "\t" + cigarLength + "\t" + alignmentLength);
                        }
						alignPos++;
						readPos++;
					}
					break;
				case SystemConfiguration.BAM_CHARD_CLIP:
				case SystemConfiguration.BAM_CPAD:
					break;
				default:
					throw new UserException("Unsupported cigar operator: " + cigarOp);
			}
		}
		return alignment;
	}

	public String getEventBases() {
		return eventBases;
	}

	public int getEventLength() {
		return eventLength;
	}
	
	/**
	 * @return the readInfo
	 */
	public AlignmentsBasic getReadInfo() {
		return readInfo;
	}

	/**
	 * @param readInfo the readInfo to set
	 */
	public void setReadInfo(AlignmentsBasic readInfo) {
		this.readInfo = readInfo;
	}

	/**
	 * @return the qpos
	 */
	public int getQpos() {
		return qpos;
	}

	/**
	 * return base in char
	 * @return
	 */
	public char getBase() {
		if(qpos < 0 || isDeletionBase())
			return 0;
		return (char) readInfo.getReadBase(qpos);
	}

	/**
	 * return base in byte
	 * @return
	 */
	public byte getByteBase() {
		return (byte)getBase();
	}

	/**
	 * return binary 4bits base
	 * @return
	 */
	public byte getBinaryBase() {
		if(qpos < 0)
			return -1;
		return readInfo.getBinaryBase(qpos);
	}

	/**
	 * return base quality
	 * @return
	 */
	public byte getBaseQuality() {
		if(qpos < 0)
			return 0;
		return readInfo.getBaseQuality(qpos);
	}

	/**
	 * return mapping quality
	 * @return
	 */
	public int getMappingQuality() {
		return readInfo.getMappingQual();
	}

	public int getPosition() {
		return readInfo.getPosition();
	}

	/**
	 * get end of alignment
	 * @return
	 */
	public int getEnd() {
		return end;
	}
	
	public int getAlignmentEnd(){
		return end+1;
	}

	public String getSample() {
		return sample;
	}

	public boolean isDeletionBase() {
		return cigarState.isDeletionBase();
	}

	public boolean isNextDeletionBase() {
		return cigarState.isNextDeletionBase();
	}

	public boolean isNextInsertBase() {
		return cigarState.isNextInsertionBase();
	}

	public boolean isNextMatchBase() {
		return cigarState.isNextMatchBase();
	}

	public boolean isInsertionAtBeginningOfRead() {
		return cigarState.isInsertionAtBeginningOfRead();
	}
}

