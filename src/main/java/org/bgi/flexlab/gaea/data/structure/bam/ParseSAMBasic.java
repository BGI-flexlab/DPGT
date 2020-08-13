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

import org.bgi.flexlab.gaea.util.CigarState;
import org.bgi.flexlab.gaea.util.SystemConfiguration;

import java.util.ArrayList;

public class ParseSAMBasic {
	
	public static String[] splitSAM(String samLine) {
		String[] alignmentArray = samLine.split("\t");
		if(alignmentArray.length < 11) {
			throw new RuntimeException("less than 11 col represent.");
		}
		return alignmentArray;
	}
	
	public static int parseFlag(String[] alignmentArray) {
		return Integer.parseInt(alignmentArray[1]);
	}
	
	public static String parseReadName(String[] alignmentArray, int flag, boolean withID) {
		String readName = "";
		if(!withID) {
			readName = alignmentArray[0];
		} else {
			if ((flag & (0x1 << 6)) != 0) {
				readName = alignmentArray[0] + "/1";
			} else if ((flag & (0x1 << 7)) != 0) {
				readName = alignmentArray[0] + "/2";
			}
		}
		return readName;
	}
	
	public static String parseChrName(String[] alignmentArray) {
		return alignmentArray[2];
	}
	
	public static int parsePosition(String[] alignmentArray, boolean is0base) {
		int position = Integer.parseInt(alignmentArray[3]);
		if(is0base) {
			position -= 1;
		}
		return position;
	}
	
	public static short parseMappingQuality(String[] alignmentArray) {
		return Short.parseShort(alignmentArray[4]);
	}
	
	//FIXME::data and method about cigar should be warraped as Class
	public static String parseCigarString(String[] alignmentArray) {
		return alignmentArray[5];
	}
	
	public static int getStartSoftClipLength(ArrayList<Integer> cigar) {
		int[] cValue = getCigarValue(cigar.get(0));
		if(cValue[0] == SystemConfiguration.BAM_CSOFT_CLIP) {
			return cValue[1];
		}
		return 0;
	}
	
	public static int getEndSoftClipLength(ArrayList<Integer> cigar) {
		int[] cValue = getCigarValue(cigar.get(cigar.size() - 1));
		if(cValue[0] == SystemConfiguration.BAM_CSOFT_CLIP) {
			return cValue[1];
		}
		return 0;
	}
	
	public static int[] getCigarValue(int cigar) {
		int[] cValue = new int[2];
		cValue[0] = (cigar & 0xf);
		cValue[1] = (cigar >> 4);
		return cValue;
	}
	
	public static int parseInsertSize(String[] alignmentArray) {
		return Integer.parseInt(alignmentArray[8]);
	}
	
	public static String parseSeq(String[] alignmentArray, int softClipStart, int softClipEnd, boolean withSoftClip) {
		if(withSoftClip) {
			return alignmentArray[9];
		} else {
			return alignmentArray[9].substring(softClipStart, alignmentArray[9].length() - softClipEnd);
		}
	}
	
	public static String parseQual(String[] alignmentArray, int softClipStart, int softClipEnd, boolean withSoftClip) {
		if(withSoftClip) {
			return alignmentArray[10];
		} else {
			return alignmentArray[10].substring(softClipStart, alignmentArray[9].length() - softClipEnd);
		}
	}
	
	public static int parseBestHitCount(String[] alignmentArray) {
		int bestHitCount = 1;
		for (int i = 11; i < alignmentArray.length; i++) {			
			if (alignmentArray[i].contains("H0:i:")) {
				bestHitCount = Integer.parseInt(alignmentArray[i].substring(5));
				break;
			} else if (alignmentArray[i].contains("X0:i:")) {
				bestHitCount = Integer.parseInt(alignmentArray[i].substring(5));
				break;
			} else if (alignmentArray[i].contains("XT:A:U")) {
				bestHitCount = 1;
			} else if (alignmentArray[i].contains("XT:A:M")) {
				bestHitCount = 2;//FIXME::not sure hit number
			} else {
				bestHitCount = 1;
			}
		}
		return bestHitCount;
	}
	
	public static String parseReadGroupID(String[] alignmentArray) {
		String rg = "";
		for (int i = 11; i < alignmentArray.length; i++) {	
			if (alignmentArray[i].contains("RG:Z:")) {
				rg = alignmentArray[i].substring(5);
				break;
			}
		}
		return rg;
	}
	
	public static int[] parseCigar(int position, CigarState cigarState) {
		int[] lenValue = new int[3];
		int end, length, baseCount;
		end = position;
		length = 0;
		baseCount = 0;
		for(int cigar : cigarState.getCigar()) {
			int[] cValue = ParseSAMBasic.getCigarValue(cigar);
			if(cValue[0] == SystemConfiguration.BAM_CMATCH) {
				end += cValue[1];
				length += cValue[1];
				baseCount += cValue[1];
			}
			if(cValue[0] == SystemConfiguration.BAM_CDEL) {
				end += cValue[1];
			}
			if(cValue[0] == SystemConfiguration.BAM_CDIFF) {
				end += cValue[1];
				length += cValue[1];
				baseCount += cValue[1];
			}
			if(cValue[0] == SystemConfiguration.BAM_CEQUAL) {
				end += cValue[1];
				length += cValue[1];
				baseCount += cValue[1];
			}
			if(cValue[0] == SystemConfiguration.BAM_CHARD_CLIP) {
				
			}
			if(cValue[0] == SystemConfiguration.BAM_CINS) {
				length += cValue[1];
			}
			if(cValue[0] == SystemConfiguration.BAM_CPAD) {
				
			}
			if(cValue[0] == SystemConfiguration.BAM_CREF_SKIP) {
				end += cValue[1];
			}
			if(cValue[0] == SystemConfiguration.BAM_CSOFT_CLIP) {
				length += cValue[1];
			}
		}
		lenValue[0] = end - 1;
		lenValue[1] = length;
		lenValue[2] = baseCount;
		
		return lenValue;
	}
}
