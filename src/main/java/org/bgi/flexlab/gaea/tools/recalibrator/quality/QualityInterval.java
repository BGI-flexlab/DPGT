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
package org.bgi.flexlab.gaea.tools.recalibrator.quality;

import org.bgi.flexlab.gaea.util.QualityUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

public class QualityInterval implements Comparable<QualityInterval> {
	public static final Set<QualityInterval> MY_EMPTY_SET = Collections.emptySet();

	public final int qStart, qEnd, fixQ, level;
	public final long obervation, error;
	public Set<QualityInterval> subIntervals;

	public QualityInterval(int qs, int qe, int fq, int lv, long ob, long errs, Set<QualityInterval> sub) {
		this.qStart = qs;
		this.qEnd = qe;
		this.fixQ = fq;
		this.level = lv;
		this.obervation = ob;
		this.error = errs;
		this.subIntervals = sub;
	}

	public QualityInterval(int qs, int lv, long ob, long errs) {
		this(qs, qs, qs, lv, ob, errs, MY_EMPTY_SET);
	}

	public QualityInterval(int qs, int qe, int lv, long ob, long errs) {
		this(qs, qe, -1, lv, ob, errs, MY_EMPTY_SET);
	}

	private boolean hasFixQuality() {
		return fixQ != -1;
	}

	public double getErrorRate() {
		if (hasFixQuality())
			return QualityUtils.qualityToErrorProbability(fixQ);
		else if (obervation == 0)
			return 0.0;
		return (error + 1) / ((double) (obervation + 1));
	}

	public byte getQual() {
		if (!hasFixQuality())
			return QualityUtils.probabilityToQuality(1 - getErrorRate(), 0);
		else
			return (byte) (fixQ & 0xff);
	}

	public boolean canMerge(QualityInterval other) {
		if ((this.qEnd + 1 != other.qStart) && (other.qEnd + 1 != this.qStart))
			return false;

		return true;
	}

	public QualityInterval merge(QualityInterval toMerge) {
		final QualityInterval left = this.compareTo(toMerge) < 0 ? this : toMerge;
		final QualityInterval right = this.compareTo(toMerge) < 0 ? toMerge : this;

		final long obs = left.obervation + right.obervation;
		final long err = left.error + right.error;

		final int level = Math.max(left.level, right.level) + 1;
		final Set<QualityInterval> subIntervals = new HashSet<QualityInterval>(Arrays.asList(left, right));
		QualityInterval merged = new QualityInterval(left.qStart, right.qEnd, -1, level, obs, err, subIntervals);

		return merged;
	}

	@Override
	public int compareTo(QualityInterval other) {
		return this.qStart - other.qStart;
	}

	public double getPenalty(int minQual) {
		return calculatePenalty(getErrorRate(), minQual);
	}

	private double calculatePenalty(final double globalErrorRate, int minQual) {
		if (globalErrorRate == 0.0)
			return 0.0;

		if (subIntervals.isEmpty()) {
			if (this.qEnd <= minQual)
				return 0;
			else {
				return (Math.abs(Math.log10(getErrorRate()) - Math.log10(globalErrorRate))) * obervation;
			}
		} else {
			double sum = 0;
			for (final QualityInterval interval : subIntervals)
				sum += interval.calculatePenalty(globalErrorRate, minQual);
			return sum;
		}
	}

	public String toString() {
		return "QQ:" + qStart + "-" + qEnd;
	}
}
