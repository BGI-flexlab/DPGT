package org.bgi.flexlab.gaea.tools.haplotypecaller.pileup;

import java.util.Collection;
import java.util.Iterator;
import java.util.NoSuchElementException;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.Utils;

public final class ReadsDownsamplingIterator implements Iterator<GaeaSamRecord>, Iterable<GaeaSamRecord> {

    private final Iterator<GaeaSamRecord> nestedReadIterator;
    private final ReadsDownsampler downsampler;
    private Iterator<GaeaSamRecord> cachedDownsampledReads = null;
    private GaeaSamRecord nextRead = null;

    /**
     * @param iter wrapped iterator from which this iterator will pull reads to be downsampled
     * @param downsampler downsampler through which the reads from the wrapped iterator will be fed
     */
    public ReadsDownsamplingIterator( Iterator<GaeaSamRecord> iter, ReadsDownsampler downsampler ) {
        Utils.nonNull(iter, "iterator must not be null");
        Utils.nonNull(downsampler, "downsampler must not be null");

        this.nestedReadIterator = iter;
        this.downsampler = downsampler;

        advanceToNextRead();
    }

    @Override
    public boolean hasNext() {
        return nextRead != null;
    }

    @Override
    public GaeaSamRecord next() {
        if ( nextRead == null ) {
            throw new NoSuchElementException("next() called when there are no more items");
        }

        final GaeaSamRecord toReturn = nextRead;
        advanceToNextRead();

        return toReturn;
    }

    private void advanceToNextRead() {
        if ( readyToReleaseReads() || fillDownsampledReadsCache() ) {
            nextRead = cachedDownsampledReads.next();
        }
        else {
            nextRead = null;
        }
    }

    private boolean readyToReleaseReads() {
        return cachedDownsampledReads != null && cachedDownsampledReads.hasNext();
    }

    private boolean fillDownsampledReadsCache() {
        while ( nestedReadIterator.hasNext() && ! downsampler.hasFinalizedItems() ) {
            downsampler.submit(nestedReadIterator.next());
        }

        if ( ! nestedReadIterator.hasNext() ) {
            downsampler.signalEndOfInput();
        }

        final Collection<GaeaSamRecord> downsampledReads = downsampler.consumeFinalizedItems();
        cachedDownsampledReads = downsampledReads.iterator();

        return cachedDownsampledReads.hasNext();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove records via a ReadsDownsamplingIterator");
    }

    @Override
    public Iterator<GaeaSamRecord> iterator() {
        return this;
    }
}
