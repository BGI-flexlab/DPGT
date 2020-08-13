package org.bgi.flexlab.gaea.tools.haplotypecaller.readfilter;

import java.io.Serializable;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;

public final class MappingQualityReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    public int minMappingQualityScore = 10;

    public Integer maxMappingQualityScore = null;

    // Command line parser requires a no-arg constructor
    public MappingQualityReadFilter() {}

    public MappingQualityReadFilter( final int minMappingQualityScore ) {
        this.minMappingQualityScore = minMappingQualityScore;
    }

    public MappingQualityReadFilter( final int minMappingQualityScore, final Integer maxMappingQualityScore) {
        this.minMappingQualityScore = minMappingQualityScore;
        this.maxMappingQualityScore = maxMappingQualityScore;
    }

    @Override
    public boolean test( final GaeaSamRecord read ) {
        final int mq = read.getMappingQuality();
        return  mq >= minMappingQualityScore
                && (maxMappingQualityScore == null || mq <= maxMappingQualityScore);
    }
}