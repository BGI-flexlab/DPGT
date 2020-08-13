package org.bgi.flexlab.gaea.tools.haplotypecaller.assembly.vertex;

import java.util.Collection;
import java.util.HashSet;

import org.bgi.flexlab.gaea.util.Utils;

final class SplitCommonSuffices extends VertexBasedTransformer {
    private final Collection<SeqVertex> alreadySplit = new HashSet<>();

    SplitCommonSuffices(final SeqGraph graph) {
        super(graph);
    }

    @Override
    boolean tryToTransform(final SeqVertex bottom) {
        Utils.nonNull(bottom);

        if (alreadySplit.contains(bottom)) {
            return false;
        } else {
            alreadySplit.add(bottom);
            return CommonSuffixSplitter.split(getGraph(), bottom);
        }
    }
}
