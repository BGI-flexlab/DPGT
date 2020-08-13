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
package org.bgi.flexlab.gaea.data.structure.bam;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Locatable;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.reads.ReadBasicCompressionInformation;
import org.bgi.flexlab.gaea.util.SystemConfiguration;

import java.io.Serializable;
import java.util.Arrays;

public class SAMCompressionInformationBasic extends ReadBasicCompressionInformation implements ParseSAMInterface, Cloneable, Serializable {
	/**
	 * flag
	 */
	protected int flag = 0;
	
	/**
	 * chromosome name index
	 */
	protected int chrNameIndex = -1;

	/**
	 * alignment start
	 */
	protected int position = 0;
	
	/**
	 * mapping quality
	 */
	protected short mappingQual = 0;
	
	/**
	 * cigar in string
	 */
	protected int[] cigars;

	public SAMCompressionInformationBasic() {
		super();
		flag = 0;
		chrNameIndex = -1;
		position = 0;
		mappingQual = 0;
	}

	public SAMCompressionInformationBasic(SAMCompressionInformationBasic read) {
		super(read);
		this.flag = read.flag;
		this.chrNameIndex = read.chrNameIndex;
		this.position = read.position;
		this.mappingQual = read.mappingQual;
		this.cigars = Arrays.copyOf(read.cigars, read.cigars.length);
	}

	public boolean parseSAM(SAMRecord samRecord) {

		flag = samRecord.getFlags();

		if(isUnmapped()) {
			return false;
		}

		chrNameIndex = samRecord.getReferenceIndex();

		position = samRecord.getAlignmentStart() - 1;

		if(position < 0) {
			return false;
		}

		mappingQual = (short) samRecord.getMappingQuality();


		if(samRecord.getCigarString().equals("*")) {
			return false;
		}

		int cigarLength = samRecord.getCigar().getCigarElements().size();
		cigars = new int[cigarLength];
		for(int i = 0; i < cigarLength; i++) {
			CigarElement cigar = samRecord.getCigar().getCigarElement(i);
			int cigarInt = (cigar.getLength() << SystemConfiguration.BAM_CIGAR_SHIFT) | CigarOperator.enumToBinary(cigar.getOperator());
			cigars[i] = cigarInt;
		}

		readBases = samRecord.getReadBases();

		qualities = samRecord.getBaseQualities();

		return true;
	}

	public boolean SAMFilter() {
		if(isUnmapped() || position < 0) {
			return false;
		}
		return true;
	}

	public int calculateReadEnd() {
		int end = position;
		for(int cigar : cigars) {
			int cigarOp = (cigar & SystemConfiguration.BAM_CIGAR_MASK);
			int cigarLength = cigar >> SystemConfiguration.BAM_CIGAR_SHIFT;
			if(cigarOp == SystemConfiguration.BAM_CMATCH || cigarOp == SystemConfiguration.BAM_CDEL ||
					cigarOp == SystemConfiguration.BAM_CREF_SKIP || cigarOp == SystemConfiguration.BAM_CEQUAL) {
				end += cigarLength;
			}
		}
		//return end;
		return end - 1;
	}

	/**
	 * get alignment end with soft clip in count
	 * @return
	 */
	public int getSoftEnd() {
		int end = calculateReadEnd();
		int cigar = cigars[cigars.length - 1];
		int cigarOp = (cigar & SystemConfiguration.BAM_CIGAR_MASK);

		if(cigarOp == SystemConfiguration.BAM_CSOFT_CLIP) {
			int cigarLength = (cigar >> SystemConfiguration.BAM_CIGAR_SHIFT);
			end += cigarLength - 1;
		}

		return end;
	}

	/**
	 * get alignment end with soft clip in count
	 * @return
	 */
	public int getSoftStart() {
		int start = position;
		int cigar = cigars[0];
		int cigarOp = (cigar & SystemConfiguration.BAM_CIGAR_MASK);

		if(cigarOp == SystemConfiguration.BAM_CSOFT_CLIP) {
			int cigarLength = cigar >> SystemConfiguration.BAM_CIGAR_SHIFT;
			start -= cigarLength;
		}

		return start;
	}


	/**
	 * hard clip read
	 * @param start length to start
	 * @param end length to end
	 */
	public void hardClip(int start, int end) {
		if(start < 0 || end < 0 || end > readBases.length - 1)
			throw new UserException("start or end < 0 or end > read length for clip read.");

		//bases
		readBases = Arrays.copyOfRange(readBases, start, end + 1);

		//qualities
		qualities = Arrays.copyOfRange(qualities, start, end + 1);
	}

	public int getNumHardClippedBases() {
		int hardClippedBasesNum = 0;
		for(int cigar : cigars) {
			int cigarOp = (cigar & SystemConfiguration.BAM_CIGAR_MASK);
			int cigarLength = cigar >> SystemConfiguration.BAM_CIGAR_SHIFT;
			if(cigarOp == SystemConfiguration.BAM_CHARD_CLIP)
				hardClippedBasesNum += cigarLength;
		}
		return hardClippedBasesNum;
	}

	public int getNumAlignedBasesCountingSoftClips(Cigar cigar) {
		int n = 0;
		if (cigar == null) return 0;

		for (final CigarElement e : cigar.getCigarElements())
			if (e.getOperator() == CigarOperator.M || e.getOperator() == CigarOperator.S)
				n += e.getLength();

		return n;
	}

	public int calcAlignmentByteArrayOffset(final Cigar cigar, final int offset, final boolean isInsertionAtBeginningOfRead, final boolean isDeletion, final int alignmentStart, final int refLocus) {
		int pileupOffset = offset;

		// Special case for reads starting with insertion
		if (isInsertionAtBeginningOfRead)
			return 0;

		// Reassign the offset if we are in the middle of a deletion because of the modified representation of the read bases
		if (isDeletion) {
			pileupOffset = refLocus - alignmentStart;
			final CigarElement ce = cigar.getCigarElement(0);
			if (ce.getOperator() == CigarOperator.S) {
				pileupOffset += ce.getLength();
			}
		}

		int pos = 0;
		int alignmentPos = 0;

		for (int iii = 0; iii < cigar.numCigarElements(); iii++) {
			final CigarElement ce = cigar.getCigarElement(iii);
			final int elementLength = ce.getLength();

			switch (ce.getOperator()) {
				case I:
				case S:
					pos += elementLength;
					if (pos >= pileupOffset) {
						return alignmentPos;
					}
					break;
				case D:
					if (!isDeletion) {
						alignmentPos += elementLength;
					} else {
						if (pos + elementLength - 1 >= pileupOffset) {
							return alignmentPos + (pileupOffset - pos);
						} else {
							pos += elementLength;
							alignmentPos += elementLength;
						}
					}
					break;
				case M:
				case EQ:
				case X:
					if (pos + elementLength - 1 >= pileupOffset) {
						return alignmentPos + (pileupOffset - pos);
					} else {
						pos += elementLength;
						alignmentPos += elementLength;
					}
					break;
				case H:
				case P:
				case N:
					break;
				default:
					throw new UserException("Unsupported cigar operator: " + ce.getOperator());
			}
		}

		return alignmentPos;
	}

	/**
	 * flag booleans
	 */
	public boolean hasMate() {
		return isQualified(SystemConfiguration.BAM_FPAIRED);
	}

	public boolean isPrpperPair() {
		return isQualified(SystemConfiguration.BAM_FPROPER_PAIR);
	}

	public boolean isUnmapped() {
		return isQualified(SystemConfiguration.BAM_FUNMAP);
	}

	public boolean isMateUnmapped() {
		return isQualified(SystemConfiguration.BAM_FMUNMAP);
	}

	public boolean isReverse() {
		return isQualified(SystemConfiguration.BAM_FREVERSE);
	}

	public boolean isMateReverse() {
		return isQualified(SystemConfiguration.BAM_FMREVERSE);
	}

	public boolean isFirstSegment() {
		return isQualified(SystemConfiguration.BAM_FREAD1);
	}

	public boolean isSecondaryAlignment() {
		return isQualified(SystemConfiguration.BAM_FSECONDARY);
	}

	public boolean isQCFailed() {
		return isQualified(SystemConfiguration.BAM_FQCFAIL);
	}

	public boolean isDup() {
		return isQualified(SystemConfiguration.BAM_FDUP);
	}

	private boolean isQualified(int config) {
		if((flag & config) != 0) {
			return true;
		}
		return false;
	}

	public boolean isInsertionAtBeginningOfRead() {
		return (cigars[0] & SystemConfiguration.BAM_CIGAR_MASK) == SystemConfiguration.BAM_CINS;
	}

	/**
	 * @return the flag
	 */
	public int getFlag() {
		return flag;
	}

	/**
	 * set flag value
	 * @param flag
	 */
	public void setFlag(int flag) {
		this.flag = flag;
	}

	/**
	 * 获取Read所在染色体的名称
	 * @return String
	 */
	public int getChrNameIndex() {
		return chrNameIndex;
	}


	/**
	 * set chromosome name index
	 * @param chrNameIndex
	 */
	public void setChrNameIndex(int chrNameIndex) {
		this.chrNameIndex = chrNameIndex;
	}

	/**
	 * 获取参考基因组上第一个碱基对的位置
	 * @return long
	 */
	public int getPosition() {
		return position;
	}


	/**
	 * set position
	 * @param position
	 */
	public void setPosition(int position) {
		this.position = position;
	}

	/**
	 * @return the mappingQual
	 */
	public short getMappingQual() {
		return mappingQual;
	}


	/**
	 * set mapping quality value
	 * @param mappingQual
	 */
	public void setMappingQual(short mappingQual) {
		this.mappingQual = mappingQual;
	}

	/**
	 * get cigars array
	 * @return
	 */
	public int[] getCigars() {
		return cigars;
	}

	/**
	 *
	 * @return
	 */
	public Cigar getCigar() {
		Cigar cigar = new Cigar();
		for(int cigarInt : cigars) {
			CigarOperator operator = CigarOperator.binaryToEnum(cigarInt & 0x0f);
			int length = cigarInt >> 4;
			CigarElement cigarElement = new CigarElement(length, operator);
			cigar.add(cigarElement);
		}

		return cigar;
	}

	public void checkCigar() {
		for (int cigar : cigars) {
			//int cigarOp = (cigar & SystemConfiguration.BAM_CIGAR_MASK);
			int cigarLength = cigar >> SystemConfiguration.BAM_CIGAR_SHIFT;
			if(cigarLength < 0) {
				throw new RuntimeException("bad cigar:" + toString());
			}
		}
	}

	/**
	 * set cigars
	 * @param cigars
	 */
	public void setCigars(int[] cigars) {
		this.cigars = cigars;
	}

	/**
	 * get cigar length
	 * @return
	 */
	public int getCigarsLength() {
		return cigars.length;
	}

	@Override
	public boolean parseSam(String samRecord) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(flag);
		sb.append("\t");
		sb.append(chrNameIndex);
		sb.append("\t");
		sb.append(position);
		sb.append("\t");
		sb.append(mappingQual);
		sb.append("\t");
		sb.append(CigarToString(cigars));
		sb.append("\t");
		sb.append(new String(readBases));
		sb.append("\t");
		sb.append(qualitiesToString(qualities));

		return sb.toString();
	}

	private String qualitiesToString(byte[] qualities) {
		byte[] newQual = new byte[qualities.length];
		for(int i = 0; i < qualities.length; i++) {
			newQual[i] = (byte) (qualities[i] + 64);
		}

		return new String(newQual);
	}


	public static String CigarToString(int[] cigars) {
		StringBuilder sb = new StringBuilder();
		for(int cigar : cigars) {
			int cigarOp = (cigar & SystemConfiguration.BAM_CIGAR_MASK);
			int cigarLength = cigar >> SystemConfiguration.BAM_CIGAR_SHIFT;
			sb.append(cigarLength);
			sb.append(SystemConfiguration.cigar2String.get(cigarOp));
		}
		return sb.toString();
	}

}
