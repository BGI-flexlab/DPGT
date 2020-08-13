package org.bgi.flexlab.gaea.tools.haplotypecaller.assembly.vertex;

import org.bgi.flexlab.gaea.util.Utils;

abstract class VertexBasedTransformer {
    /**
     * For testing purposes we sometimes want to test that can be transformed capabilities are working
     * without actually modifying the graph */
    private boolean dontModifyGraphEvenIfPossible = false;

    public boolean dontModifyGraphEvenIfPossible() { return dontModifyGraphEvenIfPossible; }
    public void setDontModifyGraphEvenIfPossible() {
        dontModifyGraphEvenIfPossible = true; }

    private final SeqGraph graph;

    VertexBasedTransformer(final SeqGraph graph){
        Utils.nonNull(graph);
        this.graph= graph;
    }

    SeqGraph getGraph() {
        return graph;
    }

    /**
     * Merge until the graph has no vertices that are candidates for merging
     */
    public boolean transformUntilComplete() {
        boolean didAtLeastOneTransform = false;
        boolean foundNodesToMerge = true;
        while( foundNodesToMerge ) {
            foundNodesToMerge = false;

            for( final SeqVertex v : graph.vertexSet() ) {
                foundNodesToMerge = tryToTransform(v);
                if ( foundNodesToMerge ) {
                    didAtLeastOneTransform = true;
                    break;
                }
            }
        }

        return didAtLeastOneTransform;
    }

    /**
     * Merge, if possible, seeded on the vertex v
     * @param v the proposed seed vertex to merge
     * @return true if some useful merging happened, false otherwise
     */
    abstract boolean tryToTransform(final SeqVertex v);
}