package org.bgi.flexlab.dpgt.utils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import java.util.ArrayList;

public class SimpleIntervalUtils {

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
        int size = Math.max((int)Math.floor((largeInterval.size() / (0.95*partitions))), 1);
        return splitIntervalBySize(largeInterval, size);
    }
}
