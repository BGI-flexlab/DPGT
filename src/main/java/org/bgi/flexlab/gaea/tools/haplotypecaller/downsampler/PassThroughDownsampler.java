package org.bgi.flexlab.gaea.tools.haplotypecaller.downsampler;

import java.util.LinkedList;
import java.util.List;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.pileup.ReadsDownsampler;
import org.bgi.flexlab.gaea.util.Utils;

public final class PassThroughDownsampler extends ReadsDownsampler {

    private LinkedList<GaeaSamRecord> selectedReads;

    public PassThroughDownsampler() {
        clearItems();
    }

    @Override
    public void submit(final GaeaSamRecord newRead ) {
        Utils.nonNull(newRead, "newRead");
        // All reads pass-through, no reads get downsampled
        selectedReads.add(newRead);
    }

    @Override
    public boolean hasFinalizedItems() {
        return ! selectedReads.isEmpty();
    }

    /**
     * Note that this list is a linked list and so doesn't support fast random access
     */
    @Override
    public List<GaeaSamRecord> consumeFinalizedItems() {
        // pass by reference rather than make a copy, for speed
        final List<GaeaSamRecord> downsampledItems = selectedReads;
        clearItems();
        return downsampledItems;
    }

    @Override
    public boolean hasPendingItems() {
        return false;
    }

    @Override
    public GaeaSamRecord peekFinalized() {
        return selectedReads.isEmpty() ? null : selectedReads.getFirst();
    }

    @Override
    public GaeaSamRecord peekPending() {
        return null;
    }

    @Override
    public int size() {
        return selectedReads.size();
    }

    @Override
    public void signalEndOfInput() {
        // NO-OP
    }

    @Override
    public void clearItems() {
        selectedReads = new LinkedList<>();
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
