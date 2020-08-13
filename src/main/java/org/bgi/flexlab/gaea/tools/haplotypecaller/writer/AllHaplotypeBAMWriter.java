package org.bgi.flexlab.gaea.tools.haplotypecaller.writer;

import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.Set;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.Haplotype;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.util.Locatable;

final class AllHaplotypeBAMWriter extends HaplotypeBAMWriter {

    public AllHaplotypeBAMWriter(final HaplotypeBAMDestination destination) {
        super(destination);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void writeReadsAlignedToHaplotypes(final Collection<Haplotype> haplotypes,
                                              final Locatable paddedReferenceLoc,
                                              final Collection<Haplotype> bestHaplotypes,
                                              final Set<Haplotype> calledHaplotypes,
                                              final ReadLikelihoods<Haplotype> readLikelihoods) {

        Utils.nonNull(haplotypes, "haplotypes cannot be null");
        Utils.nonNull(paddedReferenceLoc, "paddedReferenceLoc cannot be null");
        Utils.nonNull(readLikelihoods, "readLikelihoods cannot be null");

        writeHaplotypesAsReads(haplotypes, new LinkedHashSet<>(bestHaplotypes), paddedReferenceLoc);

        final int sampleCount = readLikelihoods.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            for (final GaeaSamRecord read : readLikelihoods.sampleReads(i)) {
                writeReadAgainstHaplotype(read);
            }
        }
    }
}
