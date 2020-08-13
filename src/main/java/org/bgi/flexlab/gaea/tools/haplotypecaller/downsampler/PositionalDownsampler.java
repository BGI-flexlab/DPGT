package org.bgi.flexlab.gaea.tools.haplotypecaller.downsampler;

import java.util.ArrayList;
import java.util.List;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.pileup.ReadsDownsampler;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.ReadCoordinateComparator;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.SAMFileHeader;

public final class PositionalDownsampler extends ReadsDownsampler {

    private final ReservoirDownsampler reservoir;

    private final SAMFileHeader header;

    private GaeaSamRecord previousRead;

    private List<GaeaSamRecord> finalizedReads;

    /**
     * Construct a PositionalDownsampler
     *
     * @param targetCoverage Maximum number of reads that may share any given alignment start position. Must be > 0
     * @param header SAMFileHeader to use to determine contig ordering. Non-null.
     */
    public PositionalDownsampler( final int targetCoverage, final SAMFileHeader header ) {
        Utils.validateArg(targetCoverage > 0, "targetCoverage must be > 0");
        Utils.nonNull(header);

        this.reservoir = new ReservoirDownsampler(targetCoverage);
        this.finalizedReads = new ArrayList<>();
        this.header = header;
        clearItems();
        resetStats();
    }

    @Override
    public void submit( final GaeaSamRecord newRead ) {
        Utils.nonNull(newRead, "newRead");

        // If we've moved to a new position, finalize the reads currently in our reservoir.
        handlePositionalChange(newRead);

        // Pass-through reads that have no assigned position, to avoid downsampling all unmapped reads
        // to the targetCoverage. Unmapped reads that do have an assigned position *will* be subject to
        // downsampling, however.
        if ( newRead.readHasNoAssignedPosition() ) {
            finalizedReads.add(newRead);
        }
        else {
            final int reservoirPreviouslyDiscardedItems = reservoir.getNumberOfDiscardedItems();
            reservoir.submit(newRead);
            incrementNumberOfDiscardedItems(reservoir.getNumberOfDiscardedItems() - reservoirPreviouslyDiscardedItems);
        }

        previousRead = newRead;
    }

    private void handlePositionalChange( final GaeaSamRecord newRead ) {
        // Use ReadCoordinateComparator to determine whether we've moved to a new start position.
        // ReadCoordinateComparator will correctly distinguish between purely unmapped reads and unmapped reads that
        // are assigned a nominal position.
        if ( previousRead != null && ReadCoordinateComparator.compareCoordinates(previousRead, newRead, header) != 0 ) {
            if ( reservoir.hasFinalizedItems() ) {
                finalizeReservoir();
            }
        }
    }

    private void finalizeReservoir() {
        finalizedReads.addAll(reservoir.consumeFinalizedItems());
        reservoir.resetStats();
    }

    @Override
    public boolean hasFinalizedItems() {
        return ! finalizedReads.isEmpty();
    }

    @Override
    public List<GaeaSamRecord> consumeFinalizedItems() {
        final List<GaeaSamRecord> toReturn = finalizedReads;
        finalizedReads = new ArrayList<>();
        return toReturn;
    }

    @Override
    public boolean hasPendingItems() {
        // The finalized items in the ReservoirDownsampler are pending items from the perspective of the
        // enclosing PositionalDownsampler
        return reservoir.hasFinalizedItems();
    }

    @Override
    public GaeaSamRecord peekFinalized() {
        return finalizedReads.isEmpty() ? null : finalizedReads.get(0);
    }

    @Override
    public GaeaSamRecord peekPending() {
        // The finalized items in the ReservoirDownsampler are pending items from the perspective of the
        // enclosing PositionalDownsampler
        return reservoir.peekFinalized();
    }

    @Override
    public int size() {
        return finalizedReads.size() + reservoir.size();
    }

    @Override
    public void signalEndOfInput() {
        finalizeReservoir();
    }

    @Override
    public void clearItems() {
        reservoir.clearItems();
        reservoir.resetStats();
        finalizedReads.clear();
        previousRead = null;
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return true;
    }

    @Override
    public void signalNoMoreReadsBefore( final GaeaSamRecord read ) {
        handlePositionalChange(read);
    }
}

