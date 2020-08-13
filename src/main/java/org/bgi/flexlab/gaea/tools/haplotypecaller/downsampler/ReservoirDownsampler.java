package org.bgi.flexlab.gaea.tools.haplotypecaller.downsampler;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.pileup.ReadsDownsampler;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.util.Utils;

public final class ReservoirDownsampler extends ReadsDownsampler {

    /**
     * size of our reservoir -- ie., the maximum number of reads from the stream that will be retained
     * (not including any undiscardable items)
     */
    private final int targetSampleSize;

    /**
     * if true, this downsampler will be optimized for the case
     * where most of the time we won't fill up anything like the
     * targetSampleSize elements.  If this is false, we will allocate
     * internal buffers to targetSampleSize initially, which minimizes
     * the cost of allocation if we often use targetSampleSize or more
     * elements.
     */
    private final boolean expectFewOverflows;

    /**
     * At times this can be a linked list or an array list, depending on how we're accessing the
     * data and whether or not we're expecting few overflows
     */
    private List<GaeaSamRecord> reservoir;

    /**
     * Are we currently using a linked list for the reservoir?
     */
    private boolean isLinkedList;

    /**
     * Count of the number of reads seen. Used by the reservoir downsampling
     * algorithm to ensure that all reads have an equal chance of making it into the reservoir.
     */
    private int totalReadsSeen;


    /**
     * Construct a ReservoirDownsampler
     *
     * @param targetSampleSize Size of the reservoir used by this downsampler.
     *
     * @param expectFewOverflows if true, this downsampler will be optimized for the case
     *                           where most of the time we won't fill up anything like the
     *                           targetSampleSize elements.  If this is false, we will allocate
     *                           internal buffers to targetSampleSize initially, which minimizes
     *                           the cost of allocation if we often use targetSampleSize or more
     *                           elements.
     */
    public ReservoirDownsampler(final int targetSampleSize, final boolean expectFewOverflows ) {
        if ( targetSampleSize <= 0 ) {
            throw new IllegalArgumentException("Cannot do reservoir downsampling with a sample size <= 0");
        }

        this.targetSampleSize = targetSampleSize;
        this.expectFewOverflows = expectFewOverflows;
        clearItems();
        resetStats();
    }

    /**
     * Construct a ReservoirDownsampler
     *
     * @param targetSampleSize Size of the reservoir used by this downsampler. Number of items retained
     *                         after downsampling will be min(totalReads, targetSampleSize)
     */
    public ReservoirDownsampler(final int targetSampleSize ) {
        this(targetSampleSize, false);
    }

    @Override
    public void submit ( final GaeaSamRecord newRead ) {
        Utils.nonNull(newRead, "newRead");

        // Only count reads that are actually eligible for discarding for the purposes of the reservoir downsampling algorithm
        totalReadsSeen++;

        if ( totalReadsSeen <= targetSampleSize ) {
            reservoir.add(newRead);
        } else {
            if ( isLinkedList ) {
                reservoir = new ArrayList<>(reservoir);
                isLinkedList = false;
            }

            final int randomSlot = GvcfMathUtils.getRandomGenerator().nextInt(totalReadsSeen);
            if ( randomSlot < targetSampleSize ) {
                reservoir.set(randomSlot, newRead);
            }
            incrementNumberOfDiscardedItems(1);
        }
    }

    @Override
    public boolean hasFinalizedItems() {
        return ! reservoir.isEmpty();
    }

    @Override
    public List<GaeaSamRecord> consumeFinalizedItems() {
        if (hasFinalizedItems()) {
            // pass reservoir by reference rather than make a copy, for speed
            final List<GaeaSamRecord> downsampledItems = reservoir;
            clearItems();
            return downsampledItems;
        } else {
            // if there's nothing here, don't bother allocating a new list
            return Collections.emptyList();
        }
    }

    @Override
    public boolean hasPendingItems() {
        return false;
    }

    @Override
    public GaeaSamRecord peekFinalized() {
        return reservoir.isEmpty() ? null : reservoir.get(0);
    }

    @Override
    public GaeaSamRecord peekPending() {
        return null;
    }

    @Override
    public int size() {
        return reservoir.size();
    }

    @Override
    public void signalEndOfInput() {
        // NO-OP
    }

    /**
     * Clear the data structures used to hold information
     */
    @Override
    public void clearItems() {
        // if we aren't expecting many overflows, allocate a linked list not an arraylist
        reservoir = expectFewOverflows ? new LinkedList<>() : new ArrayList<>(targetSampleSize);

        // it's a linked list if we allocate one
        isLinkedList = expectFewOverflows;

        // an internal stat used by the downsampling process, so not cleared by resetStats() below
        totalReadsSeen = 0;
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return false;
    }

    @Override
    public void signalNoMoreReadsBefore(final GaeaSamRecord read ) {
        Utils.nonNull(read);
    }
}

