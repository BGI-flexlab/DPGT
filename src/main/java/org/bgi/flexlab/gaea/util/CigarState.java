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
package org.bgi.flexlab.gaea.util;

import org.bgi.flexlab.gaea.data.exception.UserException;

import java.util.ArrayList;

public class CigarState {

	/**
	 * cigar操作数组，一个int表示一个cigar操作及其值，前4位为操作后28位为值 0 
	 * M alignment math(can be a sequence match or mismatch) 
	 * 1 I insertion to the reference 
	 * 2 D deletion from the reference 
	 * 3 N skipped region from the reference 
	 * 4 S soft clipping(clipped sequences present in SEQ) 
	 * 5 H hard clipping(clipped sequences NOT present in SEQ) 
	 * 6 P padding(silent deletion from padded reference) 
	 * 7 = sequence match 
	 * 8 X sequence mismatch
	 */
	private ArrayList<Integer> cigar = new ArrayList<Integer>();

	/**
	 * cigar state 
	 * 0->k: the index of the CIGAR operator that has just been processed. 
	 * 1->x: the reference coordinate of the start of k. 
	 * 2->y: the query coordinate of the start of k.
	 */
	private int[] cigarState = new int[3];

	/**
	 * flag represents the stats of the precess position. we cloud put it in the
	 * class that use this class, but it will need more functions to get the
	 * stat. 0x01 is deletion base 0x02 is next deletion base 0x04 is next
	 * insertion base
	 */
	private byte flag;

	public CigarState() {
		cigarState[0] = -1;
		cigarState[2] = 0;
		flag = 0;
	}

	public void parseCigar(String cigarv) {
		cigar = new ArrayList<Integer>();
		int begin = 0;
		for (int endv = 0; endv < cigarv.length(); endv++) {
			/*
			 * OPTION Description 
			 * M alignment math(can be a sequence match or mismatch) 
			 * I insertion to the reference D deletion from the reference 
			 * N skipped region from the reference S softclipping(clipped sequences present in SEQ)
			 * H hard clipping(clipped sequences NOT present in SEQ) 
			 * P padding(silent deletion from padded reference) 
			 * = sequence match 
			 * X sequence mismatch
			 */
			switch (cigarv.charAt(endv)) {
			case 'M': {
				begin = addCigarAndUpdateBegin(cigarv, endv, begin, 0);
				break;
			}
			case 'I': {
				begin = addCigarAndUpdateBegin(cigarv, endv, begin, 1);
				break;
			}
			case 'D': {
				begin = addCigarAndUpdateBegin(cigarv, endv, begin, 2);
				break;
			}
			case 'N': {
				begin = addCigarAndUpdateBegin(cigarv, endv, begin, 3);
				break;
			}
			case 'S': {
				begin = addCigarAndUpdateBegin(cigarv, endv, begin, 4);
				break;
			}
			case 'H': {
				begin = addCigarAndUpdateBegin(cigarv, endv, begin, 5);
				break;
			}
			case 'P': {
				begin = addCigarAndUpdateBegin(cigarv, endv, begin, 6);
				break;
			}
			case '=': {
				begin = addCigarAndUpdateBegin(cigarv, endv, begin, 7);
				break;
			}
			case 'X': {
				begin = addCigarAndUpdateBegin(cigarv, endv, begin, 8);
				break;
			}
			default:
				break;
			}
		}
	}

	public void parseCigar(int[] cigarValues) {
		for (int cigarValue : cigarValues) {
			cigar.add(cigarValue);
		}
	}

	private int addCigarAndUpdateBegin(String cigarv, int endv, int begin, int opIdx) {
		int len = Integer.parseInt(cigarv.substring(begin, endv));
		begin = endv + 1;
		int cigarValue = (len << 4) + opIdx;
		cigar.add(cigarValue);
		return begin;
	}

	public int getReadCoordWithRefCoord(int refPosition, int alignmentStartCountSoftClip) {
		cigarState[1] = alignmentStartCountSoftClip;
		for(int i = 0; i < cigar.size(); i++) {
			int op = cigar.get(i) & 0xf;
			int len = (cigar.get(i) >> 4);
			cigarState[0] = i;

			if(op == SystemConfiguration.BAM_CSOFT_CLIP || op == SystemConfiguration.BAM_CMATCH ||
					op == SystemConfiguration.BAM_CDEL || op == SystemConfiguration.BAM_CREF_SKIP ||
					op == SystemConfiguration.BAM_CEQUAL || op == SystemConfiguration.BAM_CDIFF) {
				if(refPosition >= cigarState[1] && refPosition < cigarState[1] + len) { //current cigar
					if(op == SystemConfiguration.BAM_CDEL || op == SystemConfiguration.BAM_CREF_SKIP) {
						return cigarState[2];
					} else {
						return cigarState[2] + refPosition - cigarState[1];
					}
				}

				cigarState[1] += len;
				if(op != SystemConfiguration.BAM_CDEL && op != SystemConfiguration.BAM_CREF_SKIP)
					cigarState[2] += len;
			} else {
				if(op == SystemConfiguration.BAM_CINS) {
					cigarState[2] += len;
				}
			}
		}

		throw new UserException("ref position is beyond alignment end.");
	}

	/**
	 * 解析cigar信息得到相应位点的pos
	 * 
	 * @param refPosition
	 */
	public int resolveCigar(int refPosition, int alignmentStart) {
		// determine the current CIGAR operation
		if (cigarState[0] == -1) {
			setCigarStateForNotProcessedRead(refPosition, alignmentStart);
		} else {
			int l = (cigar.get(cigarState[0]) >> 4);
			if (refPosition - cigarState[1] >= l) { // jump to the next operation
				assert (cigarState[0] < cigar.size()); // otherwise a bug: this function should not
														// be called in this case
				if (cigarState[0] + 1 >= cigar.size()) {
					throw new RuntimeException("\tcigar type:" + (cigar.get(cigarState[0]) & 0xf) + "\tcigar size:" + l
							+ "\tpos:" + refPosition + "\tpostion:" + alignmentStart);
				}

				setCigarStateForProcessedRead(l);

				assert (cigarState[0] < cigar.size()); // otherwise a bug
			} // else, do nothing
		}

		int qpos = collectPileUpInfo(refPosition);
		return qpos;
	}

	private int collectPileUpInfo(int refPosition) {
		int qpos = 0;
		if (cigarState[0] >= cigar.size()) {
			throw new RuntimeException("\tk:" + cigarState[0] + "\tx:" + cigarState[1] + "\ty:" + cigarState[2]);
		}
		int op, l;
		op = cigar.get(cigarState[0]) & 0xf;
		l = (cigar.get(cigarState[0]) >> 4);

		flag = 0; // important

		if (op == SystemConfiguration.BAM_CMATCH || op == SystemConfiguration.BAM_CEQUAL
				|| op == SystemConfiguration.BAM_CDIFF) {
			qpos = (int) (cigarState[2] + (refPosition - cigarState[1]));
		} else if (op == SystemConfiguration.BAM_CDEL || op == SystemConfiguration.BAM_CREF_SKIP) {
			qpos = -1; // FIXME: distinguish D and N!!!!!
			flag = (byte) (flag | 0x01);
		} // cannot be other operations; otherwise a bug

		if (refPosition - cigarState[1] == l - 1 && cigarState[0] + 1 < cigar.size()) {
			int nextOP = cigar.get(cigarState[0] + 1) & 0xf;
			if (nextOP == SystemConfiguration.BAM_CDEL)
				flag = (byte) (flag | 0x02);
			if (nextOP == SystemConfiguration.BAM_CINS)
				flag = (byte) (flag | 0x04);
			if (nextOP == SystemConfiguration.BAM_CMATCH || nextOP == SystemConfiguration.BAM_CEQUAL
					|| nextOP == SystemConfiguration.BAM_CDIFF) {
				flag = (byte) (flag | 0x08);
			}
		}

		return qpos;
	}

	private void setCigarStateForProcessedRead(int operatorLength) {
		int k = 0;
		int op = cigar.get(cigarState[0] + 1) & 0xf;
		if (op == SystemConfiguration.BAM_CMATCH || op == SystemConfiguration.BAM_CDEL
				|| op == SystemConfiguration.BAM_CREF_SKIP || op == SystemConfiguration.BAM_CEQUAL
				|| op == SystemConfiguration.BAM_CDIFF) { // jump to the next without a loop
			setCigarState(operatorLength);
		} else { // find the next M/D/N/=/X
			setCigarState(k, operatorLength);
		}
	}

	private void setCigarState(int cigarIndex, int length) {
		int op = cigar.get(cigarState[0]) & 0xf;
		if (op == SystemConfiguration.BAM_CMATCH || op == SystemConfiguration.BAM_CEQUAL
				|| op == SystemConfiguration.BAM_CDIFF) {
			cigarState[2] += length;
		}
		cigarState[1] += length;
		for (cigarIndex = cigarState[0] + 1; cigarIndex < cigar.size(); ++cigarIndex) {
			op = cigar.get(cigarIndex) & 0xf;
			length = cigar.get(cigarIndex) >> 4;
			if (op == SystemConfiguration.BAM_CMATCH || op == SystemConfiguration.BAM_CDEL
					|| op == SystemConfiguration.BAM_CREF_SKIP || op == SystemConfiguration.BAM_CEQUAL
					|| op == SystemConfiguration.BAM_CDIFF) {
				break;
			} else if (op == SystemConfiguration.BAM_CINS) {
				cigarState[2] += length;
			}
		}
		cigarState[0] = cigarIndex;
	}

	private void setCigarState(int l) {
		int op2 = cigar.get(cigarState[0]) & 0xf;
		if (op2 == SystemConfiguration.BAM_CMATCH || op2 == SystemConfiguration.BAM_CEQUAL
				|| op2 == SystemConfiguration.BAM_CDIFF) {
			cigarState[2] += l;
		}
		cigarState[1] += l;
		++cigarState[0];
	}

	private void setCigarStateForNotProcessedRead(int refPosition, int alignmentStart) {
		int k = 0;
		if (cigar.size() == 1) {
			setCigarStateForOneOp(alignmentStart);
		} else {
			k = findMatchOrDeletion(refPosition, alignmentStart, k);
			assert (k < cigar.size());
			cigarState[0] = k;
		}
	}

	private void setCigarStateForOneOp(int alignmentStart) {
		int cigarOp = (cigar.get(0) & 0xf);
		if (cigarOp == SystemConfiguration.BAM_CMATCH || cigarOp == SystemConfiguration.BAM_CEQUAL
				|| cigarOp == SystemConfiguration.BAM_CDIFF) {
			cigarState[0] = 0;
			cigarState[1] = alignmentStart;// 需要修改long
			cigarState[2] = 0;
		}
	}

	private int findMatchOrDeletion(int refPosition, int alignmentStart, int cigarOperatorIndex) {
		for (cigarOperatorIndex = 0, cigarState[1] = alignmentStart, cigarState[2] = 0; cigarOperatorIndex < cigar
				.size(); ++cigarOperatorIndex) {
			int op = (cigar.get(cigarOperatorIndex) & 0xf);
			int l = (cigar.get(cigarOperatorIndex) >> 4);
			if (op == SystemConfiguration.BAM_CMATCH || op == SystemConfiguration.BAM_CDEL
					|| op == SystemConfiguration.BAM_CEQUAL || op == SystemConfiguration.BAM_CDIFF) {
				if (refPosition - cigarState[1] < l) {
					break;
				} else {
					cigarState[1] += l;
					if (op == SystemConfiguration.BAM_CMATCH || op == SystemConfiguration.BAM_CEQUAL
							|| op == SystemConfiguration.BAM_CDIFF) {
						cigarState[2] += l;
					}
				}
			} else if (op == SystemConfiguration.BAM_CREF_SKIP) {
				cigarState[1] += l;
			} else if (op == SystemConfiguration.BAM_CINS || op == SystemConfiguration.BAM_CSOFT_CLIP) {
				cigarState[2] += l;
			}
		}

		return cigarOperatorIndex;
	}

	public void reSetCigarstate() {
		cigarState[0] = -1;
		cigarState[2] = 0;
	}

	public int getCurrentCigar() {
		return cigar.get(cigarState[0]);
	}

	public int getNextCigar() {
		if (cigarState[0] + 1 < cigar.size())
			return cigar.get(cigarState[0] + 1);
		else
			return -1;
	}

	public boolean isDeletionBase() {
		if ((flag & 0x01) != 0)
			return true;
		return false;
	}

	public boolean isNextDeletionBase() {
		if ((flag & 0x02) != 0)
			return true;
		return false;
	}

	public boolean isNextInsertionBase() {
		if ((flag & 0x04) != 0)
			return true;
		return false;
	}

	public boolean isNextMatchBase() {
		if ((flag & 0x08) != 0)
			return true;
		return false;
	}

	public boolean isInsertionAtBeginningOfRead() {
		int op = cigar.get(0) & 0xf;
		return op == SystemConfiguration.BAM_CINS;
	}

	public int getEventLength() {
		if(isNextDeletionBase() || isNextInsertionBase())
			return (cigar.get(cigarState[0] + 1) >> 4);
		else
			return -1;
	}

	public ArrayList<Integer> getCigar() {
		return cigar;
	}

	public int[] getCigarState() {
		return cigarState;
	}

	public void setCigar(ArrayList<Integer> cigar) {
		this.cigar = cigar;
	}

	public void setCigarState(int[] cigarState) {
		this.cigarState = cigarState;
	}
}
