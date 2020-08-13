package org.bgi.flexlab.gaea.tools.haplotypecaller.assembly.vertex;

import java.util.Arrays;

public final class SeqVertex extends BaseVertex {
    /**
     * Create a new SeqVertex with sequence and the next available id
     * @param sequence our base sequence
     */
    public SeqVertex(final byte[] sequence) {
        super(sequence);
    }

    /**
     * Create a new SeqVertex having bases of sequence.getBytes()
     * @param sequence the string representation of our bases
     */
    public SeqVertex(final String sequence) {
        super(sequence);
    }

    /**
     * Get the unique ID for this SeqVertex
     * @return a positive integer >= 0
     */
    public int getId() {
        return hashCode();
    }

    @Override
    public String toString() {
        return "SeqVertex_id_" + hashCode() + "_seq_" + getSequenceString();
    }

    /**
     * Two SeqVertex are equal only if their ids are equal
     * @param o
     * @return
     */
    @Override
    public boolean equals(final Object o) { return o == this; }

    @Override
    public int hashCode() {
        return System.identityHashCode(this);
    }

    /**
     * Return a new SeqVertex derived from this one but not including the suffix bases
     *
     * @param suffix the suffix bases to remove from this vertex
     * @return a newly allocated SeqVertex with appropriate prefix, or null if suffix removes all bases from this node
     */
    public SeqVertex withoutSuffix(final byte[] suffix) {
        final int prefixSize = sequence.length - suffix.length;
        return prefixSize > 0 ? new SeqVertex(Arrays.copyOf(sequence, prefixSize)) : null;
    }

    /**
     * Return a new SeqVertex derived from this one but not including prefix or suffix bases
     *
     * @param prefix the previx bases to remove
     * @param suffix the suffix bases to remove from this vertex
     * @return a newly allocated SeqVertex
     */
    public SeqVertex withoutPrefixAndSuffix(final byte[] prefix, final byte[] suffix) {
        final int start = prefix.length;
        final int length = sequence.length - suffix.length - prefix.length;
        final int stop = start + length;
        return length > 0 ? new SeqVertex(Arrays.copyOfRange(sequence, start, stop)) : null;
    }
}
