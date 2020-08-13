package org.bgi.flexlab.gaea.tools.haplotypecaller.assembly;

import java.util.Collections;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.ParamUtils;
import org.bgi.flexlab.gaea.util.Utils;

final class DeadEndKBestSubHaplotypeFinder implements KBestSubHaplotypeFinder {

    /**
     * Sole instance of this class.
     */
    public static final DeadEndKBestSubHaplotypeFinder INSTANCE = new DeadEndKBestSubHaplotypeFinder();

    /**
     * Prevents instantiation of more than one instance; please use {@link #INSTANCE}.
     */
    private DeadEndKBestSubHaplotypeFinder() {
    }

    @Override
    public String id() {
        return "<DEAD>";
    }

    @Override
    public String label() {
        return "&lt;DEAD&gt;";
    }

    @Override
    public Set<Pair<? extends KBestSubHaplotypeFinder, String>> subFinderLabels() {
        return Collections.emptySet();
    }

    @Override
    public int getCount() {
        return 0;
    }

    @Override
    public KBestHaplotype getKBest(final int k) {
        ParamUtils.isPositiveOrZero(k, "k cannot be negative");
        throw new IllegalArgumentException("k cannot be equal or greater to the haplotype count");
    }

    @Override
    public boolean isReference() {
        return false;
    }

    @Override
    public double score(final byte[] bases, final int offset, final int length) {
        Utils.nonNull(bases, "bases cannot be null");
        ParamUtils.isPositiveOrZero(offset, "the offset cannot be negative");
        ParamUtils.isPositiveOrZero(length, "the length cannot be negative");
        Utils.validateArg(offset + length <= bases.length, "the offset and length go beyond the array size");
        return Double.NaN;
    }
}