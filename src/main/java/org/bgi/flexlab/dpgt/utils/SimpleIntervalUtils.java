/**
This file is part of DPGT.
Copyright (C) 2022 BGI.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// License End
package org.bgi.flexlab.dpgt.utils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class SimpleIntervalUtils {
    private static final Logger logger = LoggerFactory.getLogger(SimpleIntervalUtils.class);

    /**
     * make windows of specified size by splitting large interval
     * @param largeInterval large interval(window) to be splitted
     * @param size small interval(window) size
     * @return small windows
     */
    public static ArrayList<SimpleInterval> splitIntervalBySize(final SimpleInterval largeInterval, int size) {
        ArrayList<SimpleInterval> results = new ArrayList<>();
        if (largeInterval.size() <= size) {
            results.add(largeInterval);
            return results;
        }
        int n = largeInterval.size() / size;
        for (int i = 0; i < n; ++i) {
            int s = largeInterval.getStart() + i * size;
            int e = s + size - 1;
            SimpleInterval window = new SimpleInterval(largeInterval.getContig(), s, e);
            results.add(window);
        }
        int lastEnd = results.get(n - 1).getEnd();
        if (lastEnd < largeInterval.getEnd()) {
            SimpleInterval window = new SimpleInterval(largeInterval.getContig(), lastEnd+1, largeInterval.getEnd());
            results.add(window);
        }
        return results;
    }

    /**
     * make windows by splitting large interval into number of partitions
    * @param largeInterval large interval(window) to be splitted
    * @param partitions number of partitions
    * @return small windows
    */
    public static ArrayList<SimpleInterval> splitIntervalByPartitions(final SimpleInterval largeInterval, int partitions) {
        if (partitions == 1) {
            ArrayList<SimpleInterval> result = new ArrayList<>();
            result.add(largeInterval);
            return result;
        }
        int size = Math.max((int)Math.floor((largeInterval.size() / (0.95*partitions))), 1);
        return splitIntervalBySize(largeInterval, size);
    }

    /**
     * make windows by splitting large interval into number of partitions, each small window has the similar number of bits
     * been setted
     * @param largeInterval large interval(window) to be splitted
     * @param partitions number of partitions
     * @param minVariantSites minimum number of variant sites of the small window.
     * @param bitSet bit set
     * @return small windows
     */
    public static ArrayList<SimpleInterval> splitIntervalByPartitionsAndBitSet(final SimpleInterval largeInterval, int partitions, int minVariantSites, final BitSet bitSet) {
        int totalCount = bitSet.cardinality();
        int count = Math.max((int)Math.floor((totalCount / (0.98*partitions))), minVariantSites);
        int n = 0;
        int s = bitSet.nextSetBit(0);
        int e = 0;
        ArrayList<SimpleInterval> results = new ArrayList<>();
        for (int i = s; i >= 0; i = bitSet.nextSetBit(i+1)) {
            ++n;
            if (n == count) {
                e = i;
                results.add(new SimpleInterval(largeInterval.getContig(),
                    s + largeInterval.getStart(), e + largeInterval.getStart()));
                s = bitSet.nextSetBit(i+1);
                n = 0;
            }
        }

        if (n > 0) {
            results.add(new SimpleInterval(largeInterval.getContig(), 
                s + largeInterval.getStart(), bitSet.length() + largeInterval.getStart()));
        }

        int o = checkAnyIntervalsOverlap(results);
        if (o > 0) {
            logger.error("{} overlapped {}", results.get(o-1), results.get(o));
            System.exit(1);
        }

        if (!checkIfIntervalsCoverBits(results, largeInterval, bitSet)) {
            logger.error("not all bits are covered by intervals, {}", largeInterval);
            System.exit(1);
        }

        return results;
    }

    /**
     * check if intervals cover all bits
     * @param intervals intervals
     * @param largeInterval large interval
     * @param bitSet bitset
     * @return true if intervals cover all bits
     */
    private static boolean checkIfIntervalsCoverBits(final List<SimpleInterval> intervals, final SimpleInterval largeInterval, final BitSet bitSet) {
        int numCovered = 0;
        for (SimpleInterval interval: intervals) {
            for (int i = interval.getStart(); i <= interval.getEnd(); ++i) {
                int j = i - largeInterval.getStart();
                if (bitSet.get(j)) ++numCovered;
            }
        }
        if (bitSet.cardinality() == numCovered) {
            return true;
        }
        return false;
    }

    /**
     * check if any intervals are overlapped, input intervals is assumed to be sorted and are on the same contig
     * @param intervals
     * @return index(>0) of the first interval overlapped with previous interval
     */
    private static int checkAnyIntervalsOverlap(final List<SimpleInterval> intervals) {
        if (intervals.isEmpty()) {
            return 0;
        }
        int e = intervals.get(0).getEnd();
        for (int i = 1; i < intervals.size(); ++i) {
            if (intervals.get(i).getStart() <= e) {
                return i;
            }
            e = intervals.get(i).getEnd();
        }
        return 0;
    }

}
