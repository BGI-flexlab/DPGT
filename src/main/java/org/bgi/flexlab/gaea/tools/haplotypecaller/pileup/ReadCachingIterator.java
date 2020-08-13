package org.bgi.flexlab.gaea.tools.haplotypecaller.pileup;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;

public class ReadCachingIterator implements Iterator<GaeaSamRecord> {
    
    private final Iterator<GaeaSamRecord> wrappedIter;
    private List<GaeaSamRecord> cache;
    private static final int INITIAL_CACHE_CAPACITY = 10000;

    /**
     * @param wrappedIter GaeaSamRecord iterator to wrap
     */
    public ReadCachingIterator(final Iterator<GaeaSamRecord> wrappedIter) {
        this.wrappedIter = wrappedIter;
        this.cache = new ArrayList<>(INITIAL_CACHE_CAPACITY);
    }

    @Override
    public boolean hasNext() {
        return wrappedIter.hasNext();
    }

    @Override
    public GaeaSamRecord next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("next() called when there are no more items");
        }

        final GaeaSamRecord nextRead = wrappedIter.next();
        cache.add(nextRead);
        return nextRead;
    }

    /**
     * @return All reads currently saved in the cache. The cache is emptied as a side effect of calling this.
     */
    public List<GaeaSamRecord> consumeCachedReads() {
        final List<GaeaSamRecord> oldCache = cache;
        cache = new ArrayList<>(INITIAL_CACHE_CAPACITY);
        return oldCache;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported");
    }
}

