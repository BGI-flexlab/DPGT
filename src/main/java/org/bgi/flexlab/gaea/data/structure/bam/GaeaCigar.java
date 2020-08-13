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

import java.nio.ByteBuffer;
import java.util.ArrayList;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class GaeaCigar {

	public static boolean isClipOperator(CigarOperator op) {
		return op == CigarOperator.S || op == CigarOperator.H
				|| op == CigarOperator.P;
	}

	public static Cigar unclipCigar(Cigar cigar) {
		ArrayList<CigarElement> elements = new ArrayList<CigarElement>(
				cigar.numCigarElements());
		for (CigarElement ce : cigar.getCigarElements()) {
			if (!isClipOperator(ce.getOperator()))
				elements.add(ce);
		}
		return new Cigar(elements);
	}

	public static Cigar reclipCigar(Cigar cigar, GaeaSamRecord read) {
		ArrayList<CigarElement> elements = new ArrayList<CigarElement>();
		int n = read.getCigar().numCigarElements();
		int i;
		for (i = 0; i < n; i++) {
			CigarOperator op = read.getCigar().getCigarElement(i).getOperator();
			if (isClipOperator(op))
				elements.add(read.getCigar().getCigarElement(i));
			else
				break;
		}

		elements.addAll(cigar.getCigarElements());

		for (++i; i < n; i++) {
			CigarOperator op = read.getCigar().getCigarElement(i).getOperator();
			if (isClipOperator(op))
				elements.add(read.getCigar().getCigarElement(i));
		}

		return new Cigar(elements);
	}

	public static int numberOfMatchCigarOperator(Cigar cigar) {
		int number = 0;

		if (cigar != null) {
			for (CigarElement element : cigar.getCigarElements()) {
				if (element.getOperator() == CigarOperator.M)
					number++;
			}
		}

		return number;
	}

	public static int firstIndexOfIndel(Cigar cigar) {
		int indexOfIndel = -1;

		for (int i = 0; i < cigar.numCigarElements(); i++) {
			CigarOperator op = cigar.getCigarElement(i).getOperator();
			if (op == CigarOperator.I || op == CigarOperator.D) {
				if (indexOfIndel != -1)
					return -1;
				indexOfIndel = i;
			}
		}
		return indexOfIndel;
	}

	public static Cigar moveCigarLeft(Cigar cigar, int indexOfIndel, int step) {
		if (cigar == null)
			return null;

		ArrayList<CigarElement> cigarElements = new ArrayList<CigarElement>(
				cigar.numCigarElements());
		for (int i = 0; i < indexOfIndel - 1; i++)
			cigarElements.add(cigar.getCigarElement(i));

		CigarElement ce = cigar.getCigarElement(indexOfIndel - 1);
		if (ce.getLength() < step)
			return null;
		cigarElements.add(new CigarElement(ce.getLength() - step, ce
				.getOperator()));
		cigarElements.add(cigar.getCigarElement(indexOfIndel));
		if (indexOfIndel + 1 < cigar.numCigarElements()) {
			ce = cigar.getCigarElement(indexOfIndel + 1);
			cigarElements.add(new CigarElement(ce.getLength() + step, ce
					.getOperator()));
		} else {
			cigarElements.add(new CigarElement(step, CigarOperator.M));
		}

		for (int i = indexOfIndel + 2; i < cigar.numCigarElements(); i++)
			cigarElements.add(cigar.getCigarElement(i));

		return new Cigar(cigarElements);
	}

	public static boolean cigarHasZeroSizeElement(Cigar cigar) {
		if (cigar != null) {
			for (CigarElement ce : cigar.getCigarElements()) {
				if (ce.getLength() == 0)
					return true;
			}
		}
		return false;
	}

	public static Cigar cleanCigar(Cigar cigar) {
		ArrayList<CigarElement> elements = new ArrayList<CigarElement>(
				cigar.numCigarElements() - 1);
		for (CigarElement ce : cigar.getCigarElements()) {
			if (ce.getLength() != 0
					&& (elements.size() != 0 || ce.getOperator() != CigarOperator.D)) {
				elements.add(ce);
			}
		}
		return new Cigar(elements);
	}
	
	public static Cigar decode(final ByteBuffer binaryCigar) {
        final Cigar ret = new Cigar();
        while (binaryCigar.hasRemaining()) {
            final int cigarette = binaryCigar.getInt();
            ret.add(binaryCigarToCigarElement(cigarette));
        }
        return ret;
    }
	
	private static CigarElement binaryCigarToCigarElement(final int cigarette) {
        final int binaryOp = cigarette & 0xf;
        final int length = cigarette >> 4;
        return new CigarElement(length, CigarOperator.binaryToEnum(binaryOp));
    }
}
