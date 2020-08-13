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
package org.bgi.flexlab.gaea.util;

import org.bgi.flexlab.gaea.data.structure.bam.ParseSAMBasic;
import org.bgi.flexlab.gaea.data.structure.bam.SAMInformationBasic;



/**
 * @author ZhangYong
 *
 */
public class SamRecordDatum extends SAMInformationBasic {
	/**
	 * Read ID
	 */
	private String id;

	/**
	 * Read比对到参考基因组上的次数
	 */
	private int bestHitCount;
	
	/**
	 * insert Size
	 */
	private int insertSize;
	
	/**
	 * 
	 */
	private boolean isrepeat = false;
	
	/**
	 * read end on ref, reads seq length in sam, reads base count on ref
	 */
	private int[] lenValue;
	
	/**
	 * rg index
	 */
	private int rgIndex;

	@Override
	public boolean parseSam(String samRecord) {
		String[] alignmentArray = ParseSAMBasic.splitSAM(samRecord);
		flag = ParseSAMBasic.parseFlag(alignmentArray);
		 if(isUnmapped())
			 return false;
		 id = ParseSAMBasic.parseReadName(alignmentArray, flag, true);
		 chrName = ParseSAMBasic.parseChrName(alignmentArray);
		 position = ParseSAMBasic.parsePosition(alignmentArray, true);
		 if(position < 0)
		 	return false;
		 mappingQual = ParseSAMBasic.parseMappingQuality(alignmentArray);
		 cigarString = ParseSAMBasic.parseCigarString(alignmentArray);
		 if(cigarString.equals("*"))
		  	return false;
		 cigarState = new CigarState();
	  	 cigarState.parseCigar(cigarString);
   	  	 int softClipStart = ParseSAMBasic.getStartSoftClipLength(cigarState.getCigar());
	  	 int softClipEnd = ParseSAMBasic.getEndSoftClipLength(cigarState.getCigar());
		 readSequence = ParseSAMBasic.parseSeq(alignmentArray, softClipStart, softClipEnd, true);
	  	 qualityString = ParseSAMBasic.parseQual(alignmentArray, softClipStart, softClipEnd, true);
	  	 bestHitCount = ParseSAMBasic.parseBestHitCount(alignmentArray);
	  	 insertSize = ParseSAMBasic.parseInsertSize(alignmentArray);
	  	 lenValue = ParseSAMBasic.parseCigar(position, cigarState);
	  	 return true;
	}

	/**
	 * 获取Read比对到参考基因组上的最佳比对次数
	 * @return int
	 */
	public int getBestHitCount() {
		return bestHitCount;
	}

	/**
	 * 判断Read比对到参考基因组上的次数是否为1
	 * @return true or false
	 */
	public boolean isUniqueAlignment() {
		return (bestHitCount == 1);
	}

	/**
	 * 获取Read ID
	 * @return
	 */
	public String getId() {
		return id;
	}

	/**
	 * @return the end
	 */
	public int getEnd() {
		return lenValue[0];
	}

	public int getLength() {
		return lenValue[1];
	}

	/**
	 * @return the insertSize
	 */
	public int getInsertSize() {
		return insertSize;
	}

	/**
	 * @return the isrepeat
	 */
	public boolean isRepeat() {
		return isrepeat;
	}

	/**
	 * @param isrepeat the isrepeat to set
	 */
	public void setIsrepeat(boolean isrepeat) {
		this.isrepeat = isrepeat;
	}

	/**
	 * @return the baseCount
	 */
	public int getBaseCount() {
		return lenValue[2];
	}

	/**
	 * @return the rgIndex
	 */
	public int getRgIndex() {
		return rgIndex;
	}

	/**
	 * @param rgIndex the rgIndex to set
	 */
	public void setRgIndex(int rgIndex) {
		this.rgIndex = rgIndex;
	}

	public boolean parseBamQC(String value) {
		// TODO Auto-generated method stub
		if (value.isEmpty()) {
			 return false;
		}

		String[] alignmentArray = value.split("\t");
		
		flag = Integer.parseInt(alignmentArray[0]);

		readSequence = alignmentArray[1];
		
		insertSize = Integer.parseInt(alignmentArray[2]);
	
		position = Integer.parseInt(alignmentArray[3]);
	
		if(position < 0 ) {
			return false;
		}
		cigarString = alignmentArray[4];
		cigarState = new CigarState();
		cigarState.parseCigar(cigarString);
		bestHitCount = Integer.parseInt(alignmentArray[5]);
		if(Integer.parseInt(alignmentArray[6]) == 1) {
			isrepeat = true;
		}
		mappingQual = Short.parseShort(alignmentArray[7]);
		if(alignmentArray.length >= 9)
		rgIndex = Integer.parseInt(alignmentArray[8]);
		if(alignmentArray.length >= 10)
			qualityString = alignmentArray[9];
		lenValue = ParseSAMBasic.parseCigar(position, cigarState);
		return true;
	}
}
