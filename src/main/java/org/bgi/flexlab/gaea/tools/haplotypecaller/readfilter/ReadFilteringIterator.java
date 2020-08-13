package org.bgi.flexlab.gaea.tools.haplotypecaller.readfilter;

import java.util.Iterator;
import java.util.NoSuchElementException;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.Utils;

public class ReadFilteringIterator implements Iterator<GaeaSamRecord>, Iterable<GaeaSamRecord> {

    private final Iterator<GaeaSamRecord> nestedIterator;
    private final ReadFilter readFilter;
    private GaeaSamRecord nextRead;

    /**
     * Create a ReadFilteringIterator given a pre-existing iterator of reads and a read filter.
     * Only reads that pass the filter will be returned from this iterator.
     *
     * @param nestedIterator underlying iterator from which to pull reads (may not be null)
     * @param readFilter filter to apply to the reads (may not be null)
     */
    public ReadFilteringIterator( final Iterator<GaeaSamRecord> nestedIterator, final ReadFilter readFilter ) {
        Utils.nonNull(nestedIterator);
        Utils.nonNull(readFilter);

        this.nestedIterator = nestedIterator;
        this.readFilter = readFilter;
        this.nextRead = loadNextRead();
    }

    @Override
    public boolean hasNext() {
        return nextRead != null;
    }

    @Override
    public GaeaSamRecord next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("Iterator exhausted");
        }

        final GaeaSamRecord toReturn = nextRead;
        nextRead = loadNextRead();
        return toReturn;
    }

    private GaeaSamRecord loadNextRead() {
        while ( nestedIterator.hasNext() ) {
            final GaeaSamRecord candidate = nestedIterator.next();
            if ( readFilter.test(candidate) ) {
                return candidate;
            }
        }
        return null;
    }

    @Override
    public Iterator<GaeaSamRecord> iterator() {
        return this;
    }
}

