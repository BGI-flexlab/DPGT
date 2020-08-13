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

import org.bgi.flexlab.gaea.data.structure.reads.ReadInformation;
import org.bgi.flexlab.gaea.util.CigarState;
import org.bgi.flexlab.gaea.util.SystemConfiguration;

public abstract class SAMInformationBasic extends ReadInformation implements ParseSAMInterface{
	/**
	 * flag
	 */
	protected int flag;
	
	/**
	 * 染色体名称
	 */
	protected String chrName;

	/**
	 * 参考基因组上第一个碱基对的坐标值
	 */
	protected int position;
	
	/**
	 * mapping quality
	 */
	protected short mappingQual = 0;
	
	/**
	 * cigar in string
	 */
	protected String cigarString;
	
	/**
	 * cigar state
	 */
	protected CigarState cigarState;
	
	@Override
	public boolean parseSam(String samline) {
		
		String[] alignmentArray = ParseSAMBasic.splitSAM(samline);
		
		flag = ParseSAMBasic.parseFlag(alignmentArray);
		
		if(isUnmapped()) {
			return false;
		}
		
		chrName = ParseSAMBasic.parseChrName(alignmentArray);
		
		position = ParseSAMBasic.parsePosition(alignmentArray, true);
		
		if(position < 0) {
			return false;
		}
		
		mappingQual = ParseSAMBasic.parseMappingQuality(alignmentArray);
		
		cigarString = ParseSAMBasic.parseCigarString(alignmentArray);
		
		if(cigarString.equals("*")) {
			return false;
		}
		
		cigarState = new CigarState();
		
		cigarState.parseCigar(cigarString);
		
		int softClipStart = ParseSAMBasic.getStartSoftClipLength(cigarState.getCigar());
		
		int softClipEnd = ParseSAMBasic.getEndSoftClipLength(cigarState.getCigar());
		
		readSequence = ParseSAMBasic.parseSeq(alignmentArray, softClipStart, softClipEnd, false);
		
		qualityString = ParseSAMBasic.parseQual(alignmentArray, softClipStart, softClipEnd, false);
				
		return true;
	}

	//@Override
	//public boolean parseBamQC(String samRecord) {
	//	return true;
	//}
	
	@Override
	public boolean SAMFilter() {
		if(isUnmapped() || cigarString.equals("*") || readSequence.length() > qualityString.length() || position < 0) {
			return false;
		}
		return true;
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

	/**
	 * 获取参考基因组上第一个碱基对的位置
	 * @return long
	 */
	public int getPosition() {
		return position;
	}

	/**
	 * 获取Read所在染色体的名称
	 * @return String
	 */
	public String getChrName() {
		return chrName;
	}

	public String getCigarString() {
		return cigarString;
	}

	/**
	 * @return the flag
	 */
	public int getFlag() {
		return flag;
	}

	/**
	 * @return the mappingQual
	 */
	public short getMappingQual() {
		return mappingQual;
	}
	
	public CigarState getCigarState() {
		return cigarState;
	}
}

