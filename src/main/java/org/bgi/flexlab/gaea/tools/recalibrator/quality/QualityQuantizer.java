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

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.util.QualityUtils;

import java.util.*;

public class QualityQuantizer {
	private List<Long> observationNumbers = null;
	private final int level;
	private final int minimumQuality;

	public QualityQuantizer(final List<Long> obser, final int level, final int minQuality) {
		this.observationNumbers = obser;
		this.level = level;
		this.minimumQuality = minQuality;

		if (Collections.min(obser) < 0)
			throw new UserException("Quality score has negative values!");
		if (level < 0)
			throw new UserException("quality level must be >= 0!");
		if (minQuality < 0)
			throw new UserException("minimum quality must be >= 0!");
	}

	private TreeSet<QualityInterval> quantize() {
		final TreeSet<QualityInterval> intervals = new TreeSet<QualityInterval>();
		for (int qStart = 0; qStart < observationNumbers.size(); qStart++) {
			final long nObs = observationNumbers.get(qStart);
			final double errorRate = QualityUtils.qualityToErrorProbability((byte) qStart);
			final double nErrors = nObs * errorRate;
			final QualityInterval qi = new QualityInterval(qStart, 0, nObs, (int) Math.floor(nErrors));
			intervals.add(qi);
		}

		while (intervals.size() > level) {
			mergeLowestPenaltyIntervals(intervals);
		}

		return intervals;
	}

	private void mergeLowestPenaltyIntervals(final TreeSet<QualityInterval> intervals) {
		final Iterator<QualityInterval> iter = intervals.iterator();
		final Iterator<QualityInterval> iterator = intervals.iterator();
		iterator.next();

		QualityInterval minMerge = null;
		while (iterator.hasNext()) {
			final QualityInterval left = iter.next();
			final QualityInterval right = iterator.next();
			if (left.canMerge(right)) {
				final QualityInterval merged = left.merge(right);
				if (minMerge == null || (merged.getPenalty(minimumQuality) < minMerge.getPenalty(minimumQuality))) {
					minMerge = merged;
				}
			}
		}

		intervals.removeAll(minMerge.subIntervals);
		intervals.add(minMerge);
	}

	public List<Byte> getIntervals() {
		TreeSet<QualityInterval> qualityIntervals = quantize();
		final List<Byte> map = new ArrayList<Byte>(observationNumbers.size());
		map.addAll(Collections.nCopies(observationNumbers.size(), Byte.MIN_VALUE));
		for (final QualityInterval interval : qualityIntervals) {
			for (int q = interval.qStart; q <= interval.qEnd; q++) {
				map.set(q, interval.getQual());
			}
		}

		if (Collections.min(map) == Byte.MIN_VALUE)
			throw new UserException("quantized quality score map contains an un-initialized value!");

		return map;
	}
}
