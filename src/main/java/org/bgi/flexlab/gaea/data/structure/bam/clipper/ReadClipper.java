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
package org.bgi.flexlab.gaea.data.structure.bam.clipper;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.bam.clipper.algorithm.ReadClippingAlgorithm;
import org.bgi.flexlab.gaea.util.CigarState;
import org.bgi.flexlab.gaea.util.SystemConfiguration;

import java.util.ArrayList;

public class ReadClipper {
	private boolean isClip;
	private final GaeaSamRecord read;
	private ArrayList<ClippingRegion> crs = null;

	public ReadClipper(GaeaSamRecord read) {
		this.read = read;
		isClip = false;
		crs = new ArrayList<ClippingRegion>();
	}

	public GaeaSamRecord getRead() {
		return this.read;
	}

	public boolean isClipper() {
		return isClip;
	}

	public GaeaSamRecord clipRead(ReadClippingAlgorithm algorithm) {
		if (crs.size() == 0)
			return read;

		isClip = true;

		GaeaSamRecord clippedRead = read;
		for (ClippingRegion op : crs) {
			final int readLength = clippedRead.getReadLength();

			if (op.getStart() < readLength) {
                ClippingRegion fixedOperation = op;
                if (op.getStop() >= readLength)
                    fixedOperation = new ClippingRegion(op.getStart(), readLength - 1);

                clippedRead = algorithm.apply(clippedRead,fixedOperation);
            }
		}
		isClip = true;
		crs.clear();
		if (clippedRead.isEmpty())
			return GaeaSamRecord.emptyRead(clippedRead);
		return clippedRead;
	}

	public GaeaSamRecord clipLowQualityEnds(byte lowQuality, ReadClippingAlgorithm algorithm) {
		if (read.isEmpty())
			return read;

		int leftIndex = 0;
		int readLength = read.getReadLength();
		int rightIndex = readLength - 1;

		byte[] qualities = read.getBaseQualities();

		while (leftIndex <= rightIndex && qualities[leftIndex] <= lowQuality)
			leftIndex++;
		while (rightIndex >= 0 && qualities[rightIndex] <= lowQuality)
			rightIndex--;

		if (rightIndex < leftIndex)
			return GaeaSamRecord.emptyRead(read);

		if(leftIndex >= readLength)
			leftIndex = readLength - 1;
		
		if(rightIndex < 0)
			rightIndex = 0;
		
		if (rightIndex != (read.getReadLength() - 1)){
			crs.add(new ClippingRegion(rightIndex, readLength - 1));
		}
		
		if (leftIndex != 0){
			crs.add(new ClippingRegion(0, leftIndex));
		}

		return clipRead(algorithm);
	}

	public static AlignmentsBasic hardClipByReferenceCoordinatesLeftTail ( AlignmentsBasic read, int start) {
		return hardClipByReferenceCoordinates(read, start, -1);
	}

	public static AlignmentsBasic hardClipByReferenceCoordinatesRightTail (AlignmentsBasic read, int end) {
		return hardClipByReferenceCoordinates(read, -1, end);
	}

	public static AlignmentsBasic hardClipLowQualEnds(AlignmentsBasic read, byte lowQual) {
		if (read.isEmpty())
			return read;

		final byte [] quals = read.getQualities();
		final int readLength = read.getReadLength();
		int leftClipIndex = 0;
		int rightClipIndex = readLength - 1;

		// check how far we can clip both sides
		while (rightClipIndex >= 0 && quals[rightClipIndex] <= lowQual) rightClipIndex--;
		while (leftClipIndex < readLength && quals[leftClipIndex] <= lowQual) leftClipIndex++;

		if(rightClipIndex == readLength - 1 && leftClipIndex == 0)
			return read;
		AlignmentsBasic clippedRead = null;
		// if the entire read should be clipped, then return an empty read.
		if (leftClipIndex > rightClipIndex) {
			clippedRead = new AlignmentsBasic(read);
			clippedRead.emptyRead();
		}

		int start = 0, end = readLength - 1;
		if (rightClipIndex < readLength - 1) {
			end = rightClipIndex;
		}
		if (leftClipIndex > 0 ) {
			start = leftClipIndex;
		}
		clippedRead = hardClipByReadCoordinates(read, start, end);
		return clippedRead;
	}

	/**
	 * Generic functionality to hard clip a read, used internally by hardClipByReferenceCoordinatesLeftTail
	 * and hardClipByReferenceCoordinatesRightTail. Should not be used directly.
	 *
	 * Note, it REQUIRES you to give the directionality of your hard clip (i.e. whether you're clipping the
	 * left of right tail) by specifying either refStart < 0 or refStop < 0.
	 *
	 * @param refStart  first base to clip (inclusive)
	 * @param refStop last base to clip (inclusive)
	 * @return a new read, without the clipped bases
	 */
	//@Requires({"!read.getReadUnmappedFlag()", "refStart < 0 || refStop < 0"})  // can't handle unmapped reads, as we're using reference coordinates to clip
	protected static AlignmentsBasic hardClipByReferenceCoordinates(AlignmentsBasic read, int refStart, int refStop) {
		if (read.isEmpty() || read.isUnmapped())
			return read;

		int start;
		int stop;
		CigarState cigarState = new CigarState();
		cigarState.parseCigar(read.getCigars());

		// Determine the read coordinate to start and stop hard clipping
		if (refStart < 0) {
			if (refStop < 0)
				throw new UserException("Only one of refStart or refStop must be < 0, not both (" + refStart + ", " + refStop + ")");
			start = 0;
			stop = cigarState.getReadCoordWithRefCoord(refStop, read.getSoftStart());
		}
		else {
			if (refStop >= 0)
				throw new UserException("Either refStart or refStop must be < 0 (" + refStart + ", " + refStop + ")");
			start = cigarState.getReadCoordWithRefCoord(refStart, read.getSoftStart());
			stop = read.getReadLength() - 1;
		}

		if (start < 0 || stop > read.getReadLength() - 1)
			throw new UserException("Trying to clip before the start or after the end of a read");

		if ( start > stop )
			throw new UserException(String.format("START (%d) > (%d) STOP -- this should never happen -- call Mauricio!", start, stop));

		if ( start > 0 && stop < read.getReadLength() - 1)
			throw new UserException(String.format("Trying to clip the middle of the read: start %d, stop %d", start, stop));

		if(start == 0 && stop == read.getReadLength() - 1)
			return read;

		//System.err.println("clipped read ref start:" + refStart + "\tclipped ref read end:" + refStop);
		//System.err.println("clipped read start:" + start + "\tclipped read end:" + stop);

		AlignmentsBasic clippedRead = new AlignmentsBasic(read);

		int[] newCigars;
		//clip start
		if(refStart > 0 ) {
			//base and quality
			clippedRead.hardClip(start, read.getReadLength() - 1);

			//position
			if(refStart > clippedRead.getPosition())
				clippedRead.setPosition(refStart);

			//cigar
			int currentCigarIndex = cigarState.getCigarState()[0];
			newCigars = new int[read.getCigars().length - currentCigarIndex];
			int currentCigar = cigarState.getCurrentCigar();

			if((currentCigar & 0xf) == SystemConfiguration.BAM_CDEL) {
				newCigars[0] = ((currentCigar & 0xf) | ((refStart - cigarState.getCigarState()[1]) << 4));
			} else
				newCigars[0] = ((currentCigar & 0xf) | ((start - cigarState.getCigarState()[2]) << 4));
			for(int i = currentCigarIndex + 1; i < read.getCigars().length; i++) {
				newCigars[i - currentCigarIndex] = read.getCigars()[i];
			}
			clippedRead.setCigars(newCigars);
		} else if(refStop > 0 ) {
			//base and quality
			clippedRead.hardClip(0, stop);
			//no need to deal with position

			//cigar
			int currentCigarIndex = cigarState.getCigarState()[0];
			newCigars = new int[currentCigarIndex + 1];
			int currentCigar = cigarState.getCurrentCigar();

			for(int i = 0; i < currentCigarIndex; i++) {
				newCigars[i] = read.getCigars()[i];
			}
			if((currentCigar & 0xf) == SystemConfiguration.BAM_CDEL) {
				newCigars[currentCigarIndex] = ((currentCigar & 0xf) | (refStop - cigarState.getCigarState()[1] << 4));
			} else
				newCigars[currentCigarIndex] = ((currentCigar & 0xf) | (stop - cigarState.getCigarState()[2] << 4));
			clippedRead.setCigars(newCigars);
		}
		//System.err.println("read:" + read.toString());
		//System.err.println("clipped read:" + clippedRead.toString());

		return clippedRead;
		//return hardClipByReadCoordinates(read, start, stop);
	}

	protected static AlignmentsBasic hardClipByReadCoordinates(AlignmentsBasic read, int readStart, int readStop) {
		if (read.isEmpty() || read.isUnmapped() || (readStart == 0 && readStop == read.getReadLength() - 1))
			return read;

		if (readStart < 0 || readStop > read.getReadLength() - 1)
			throw new UserException("Trying to clip before the start or after the end of a read");

		if ( readStart > readStop )
			throw new UserException(String.format("START (%d) > (%d) STOP -- this should never happen -- call Mauricio!", readStart, readStop));

		//if ( readStart > 0 && readStop < read.getReadLength() - 1)
		//	throw new UserException(String.format("Trying to clip the middle of the read: start %d, stop %d", readStart, readStop));

		//System.err.println("clipped read start:" + readStart + "\tclipped read end:" + readStop);

		AlignmentsBasic clippedRead = new AlignmentsBasic(read);

		//clipped bases and qualities
		clippedRead.hardClip(readStart, readStop);

		//deal with cigars
		if(read.getCigars().length == 1 && (read.getCigars()[0] & 0x0f) == SystemConfiguration.BAM_CMATCH) {//deal with the normal case
			int newCigarLen = readStop - readStart + 1;
			int[] cigars = new int[1];
			cigars[0] = ((newCigarLen << 4));
			clippedRead.setCigars(cigars);

			clippedRead.setPosition(clippedRead.getPosition() + readStart);

			//System.err.println("read:" + read.toString());
			//System.err.println("clipped read:" + clippedRead.toString());
			return clippedRead;
		}

		ArrayList<Integer> newCigars = new ArrayList<>();
		int cigarForwardIndex = -1;
		int cigarReadPosition = 0;
		int position = read.getPosition();
		for(int i = 0; i < read.getCigars().length; i++) {
			int cigar = read.getCigars()[i];
			int op = cigar & 0xf;
			int len = (cigar >> 4);

			if(op == SystemConfiguration.BAM_CSOFT_CLIP || op == SystemConfiguration.BAM_CMATCH ||
					op == SystemConfiguration.BAM_CINS  || op == SystemConfiguration.BAM_CEQUAL ||
					op == SystemConfiguration.BAM_CDIFF) {

				if(cigarReadPosition + len > readStart) { //cigar beyond read start
					if(cigarForwardIndex == -1) {// meet read start
						int newCigarLen = len - (readStart - cigarReadPosition);
						int newCigar = (op | (newCigarLen << 4));
						newCigars.add(newCigar);
						cigarForwardIndex = i;

						if(op != SystemConfiguration.BAM_CSOFT_CLIP && op != SystemConfiguration.BAM_CINS) {
							position += readStart - cigarReadPosition;
						}

						if(cigarReadPosition + len > readStop) {// same cigar as the read start
							newCigars.remove(newCigars.size() - 1);
							newCigarLen = readStop - readStart + 1;
							newCigar = (op | (newCigarLen << 4));
							newCigars.add(newCigar);
							break;
						}
					} else {
						if (cigarReadPosition + len <= readStop) //have not meet read stop
							newCigars.add(cigar);
						else { //meet read stop
							int newCigarLen = readStop - cigarReadPosition + 1;
							int newCigar = (op | (newCigarLen << 4));
							newCigars.add(newCigar);
							break;
						}
					}
				}
				cigarReadPosition += len;
			} else
				newCigars.add(cigar);

			if (cigarForwardIndex == -1 && (op == SystemConfiguration.BAM_CDEL || op == SystemConfiguration.BAM_CMATCH
					|| op == SystemConfiguration.BAM_CEQUAL || op == SystemConfiguration.BAM_CDIFF)) {
				position += len;
			}
		}
		clippedRead.setPosition(position);

		int[] finalCigars = new int[newCigars.size()];
		int index = 0;
		for(int cigar : newCigars) {
			finalCigars[index++] = cigar;
		}
		clippedRead.setCigars(finalCigars);

		//System.err.println("read:" + read.toString());
		//System.err.println("clipped read:" + clippedRead.toString());
		return clippedRead;
	}

}
