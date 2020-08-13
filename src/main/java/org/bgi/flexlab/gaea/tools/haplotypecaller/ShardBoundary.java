package org.bgi.flexlab.gaea.tools.haplotypecaller;

import java.io.Serializable;
import java.util.Objects;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.util.Locatable;

public class ShardBoundary implements Locatable, Serializable {
    private static final long serialVersionUID = 1L;

    private final GenomeLocation interval;
    private final GenomeLocation paddedInterval;

    /**
     * Create a new ShardBoundary from the given intervals
     *
     * @param interval the interval covered by the shard
     * @param paddedInterval the interval covered by the shard's padding, must contain the shard interval
     */
    public ShardBoundary(final GenomeLocation interval, final GenomeLocation paddedInterval) {
        Utils.nonNull(interval);
        Utils.nonNull(paddedInterval);
        //Utils.validateArg(paddedInterval.contains(interval), "interval must be contained within paddedInterval");
        this.interval = interval;
        this.paddedInterval = paddedInterval;
    }


    @Override
    public String getContig() {
        return interval.getContig();
    }

    /**
     * @return start of the shard boundary without padding
     */
    @Override
    public int getStart() {
        return interval.getStart();
    }

    /**
     * @return end of the shard boundary without padding
     */
    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    /**
     * @return the interval this boundary covers without padding
     */
    public GenomeLocation getInterval() {
        return interval;
    }

    /**
     * @return the interval this boundary covers including padding
     */
    public GenomeLocation getPaddedInterval() {
        return paddedInterval;
    }

    @Override
    public boolean equals(Object o) {
        if ( this == o ) {
            return true;
        }
        if ( o == null || getClass() != o.getClass() ) {
            return false;
        }

        final ShardBoundary key = (ShardBoundary) o;

        if ( !Objects.equals(interval, key.interval) ) {
            return false;
        }
        return Objects.equals(paddedInterval, key.paddedInterval);

    }

    @Override
    public int hashCode() {
        return Objects.hash(interval, paddedInterval);
    }
}
