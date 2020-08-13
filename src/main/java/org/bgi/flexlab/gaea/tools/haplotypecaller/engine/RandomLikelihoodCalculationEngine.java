package org.bgi.flexlab.gaea.tools.haplotypecaller.engine;

import java.util.List;
import java.util.Map;
import java.util.Random;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.Haplotype;
import org.bgi.flexlab.gaea.tools.haplotypecaller.LikelihoodMatrix;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.haplotypecaller.SampleList;
import org.bgi.flexlab.gaea.tools.haplotypecaller.allele.AlleleList;
import org.bgi.flexlab.gaea.tools.haplotypecaller.allele.IndexedAlleleList;
import org.bgi.flexlab.gaea.tools.haplotypecaller.assembly.AssemblyResultSet;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.util.Utils;

public final class RandomLikelihoodCalculationEngine implements ReadLikelihoodCalculationEngine {

    @Override
    public ReadLikelihoods<Haplotype> computeReadLikelihoods(final AssemblyResultSet assemblyResultSet,
                                                  final SampleList samples,
                                                  final Map<String, List<GaeaSamRecord>> reads) {
        Utils.nonNull(assemblyResultSet, "assemblyResultSet is null");
        Utils.nonNull(samples, "samples is null");
        Utils.nonNull(reads, "perSampleReadList is null");
        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(assemblyResultSet.getHaplotypeList());
        final ReadLikelihoods<Haplotype> result = new ReadLikelihoods<>(samples, haplotypes, reads);
        final Random rnd = GvcfMathUtils.getRandomGenerator();
        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = haplotypes.numberOfAlleles();
        for (int i = 0; i < sampleCount; i++)  {
            final List<GaeaSamRecord> sampleReads = result.sampleReads(i);
            final int readCount = sampleReads.size();
            final LikelihoodMatrix<Haplotype> sampleLikelihoods = result.sampleMatrix(i);
            for (int a = 0; a < alleleCount; a++) {
                for (int r = 0; r < readCount; r++) {
                    sampleLikelihoods.set(a, r, -Math.abs(rnd.nextDouble()));
                }
            }
        }
        return result;
    }

    @Override
    public void close() {
    }

}