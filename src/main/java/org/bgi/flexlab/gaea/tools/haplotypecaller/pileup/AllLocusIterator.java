package org.bgi.flexlab.gaea.tools.haplotypecaller.pileup;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.NoSuchElementException;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.util.PeekableIterator;

public class AllLocusIterator implements Iterator<AlignmentContext> {
    private final PeekableIterator<AlignmentContext> nestedLocusIterator;
    private final GenomeLocation interval;

    private int currentPosition;
    private AlignmentContext nextPileup;

    /**
     * @param interval The single interval over whose loci we'll be iterating
     * @param nestedLocusIterator Provider of AlignmentContexts that may lie within the interval. Must return AlignmentContexts
     *                            that are on the same contig as the provided interval.
     */
    public AllLocusIterator(final GenomeLocation interval, final Iterator<AlignmentContext> nestedLocusIterator) {
        Utils.nonNull(interval);
        Utils.nonNull(nestedLocusIterator);
        
        this.nestedLocusIterator = new PeekableIterator<>(nestedLocusIterator);
        this.interval = interval;
        this.currentPosition = interval.getStart();
//        this.currentPosition = -1;

        // Sanity check:
        if ( this.nestedLocusIterator.peek() != null && ! this.nestedLocusIterator.peek().getContig().equals(interval.getContig()) ) {
            throw new IllegalArgumentException("Locus iterator must be over the same contig as the interval provided");
        }

        nextPileup = advance();
    }

    @Override
    public boolean hasNext() {
        return nextPileup != null;
    }

    @Override
    public AlignmentContext next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("next() called when there are no more items");
        }

        final AlignmentContext toReturn = nextPileup;
        nextPileup = advance();
        return toReturn;
    }

    private AlignmentContext advance() {
        // If we're out of loci, pull on the nested locus iterator until it's exhausted (caller may be relying on this),
        // and then return null
        if ( currentPosition > interval.getEnd() ) {
            while ( nestedLocusIterator.hasNext() ) {
                nestedLocusIterator.next();
            }

            return null;
        }

        // Discard AlignmentContexts from the nested locus iterator that are before our current position
        while ( nestedLocusIterator.hasNext() && nestedLocusIterator.peek().getStart() < currentPosition ) {
            nestedLocusIterator.next();
        }

        final AlignmentContext nextNestedPileup = nestedLocusIterator.peek();
        AlignmentContext toReturn;
        
        if( nextNestedPileup != null ) {
            if ( nextNestedPileup.getStart() == currentPosition ) {
                toReturn = nestedLocusIterator.next();
            }else {
                toReturn = createEmptyAlignmentContextForPosition(currentPosition);
            }
        }else
            toReturn = createEmptyAlignmentContextForPosition(currentPosition);

//        // No more pileups from the nested iterator? then always return an empty pileup
//        if ( nextNestedPileup == null ) {
//            toReturn = createEmptyAlignmentContextForPosition(currentPosition);
//
//        // If the pileup from the nested iterator matches our current position, return it
//        } else if ( nextNestedPileup.getStart() == currentPosition ) {
//            toReturn = nestedLocusIterator.next();
//
//        // Otherwise, the next pileup from our nested iterator must come after our current position,
//        // so keep it around and return an empty pileup for the current position
//        } else {
//            toReturn = createEmptyAlignmentContextForPosition(currentPosition);
//        }

        currentPosition++;
        return toReturn;
    }

    private AlignmentContext createEmptyAlignmentContextForPosition(final int position) {
        final GenomeLocation positionInterval = new GenomeLocation(interval.getContig(), position, position);
        return new AlignmentContext(positionInterval, new ReadPileup(positionInterval, new ArrayList<>()));
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported");
    }
}

