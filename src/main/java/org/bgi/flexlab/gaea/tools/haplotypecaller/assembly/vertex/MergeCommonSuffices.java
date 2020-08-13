package org.bgi.flexlab.gaea.tools.haplotypecaller.assembly.vertex;

import org.bgi.flexlab.gaea.util.Utils;

final class MergeCommonSuffices extends VertexBasedTransformer {
    MergeCommonSuffices(final SeqGraph graph) {
        super(graph);
    }

    @Override
    boolean tryToTransform(final SeqVertex bottom) {
        Utils.nonNull(bottom);
        return SharedSequenceMerger.merge(getGraph(), bottom);
    }
}
