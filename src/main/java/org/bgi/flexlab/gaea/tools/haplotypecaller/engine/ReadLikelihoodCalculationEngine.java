package org.bgi.flexlab.gaea.tools.haplotypecaller.engine;

import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.Haplotype;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.haplotypecaller.SampleList;
import org.bgi.flexlab.gaea.tools.haplotypecaller.assembly.AssemblyResultSet;

public interface ReadLikelihoodCalculationEngine {

    enum Implementation {
        /**
         * Classic full pair-hmm all haplotypes vs all reads.
         */
        PairHMM,

        /**
         * Random likelihoods, used to establish a baseline benchmark for other meaningful implementations.
         */
        Random
    }


    /**
     * Calculates the likelihood of reads across many samples evaluated against haplotypes resulting from the
     * active region assembly process.
     *
     * @param assemblyResultSet the input assembly results.
     * @param samples the list of targeted samples.
     * @param perSampleReadList the input read sets stratified per sample.
     *
     * @throws IllegalArgumentException if any parameter is {@code null}.
     *
     * @return never {@code null}, and with at least one entry for input sample (keys in {@code perSampleReadList}.
     *    The value maps can be potentially empty though.
     */
    public ReadLikelihoods<Haplotype> computeReadLikelihoods(AssemblyResultSet assemblyResultSet, SampleList samples,
                                                             Map<String, List<GaeaSamRecord>> perSampleReadList);

    /**
     * This method must be called when the client is done with likelihood calculations.
     * It closes any open resources.
     */
    public void close();
}
