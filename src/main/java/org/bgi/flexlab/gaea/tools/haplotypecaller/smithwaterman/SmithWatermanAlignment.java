package org.bgi.flexlab.gaea.tools.haplotypecaller.smithwaterman;

import htsjdk.samtools.Cigar;

public interface SmithWatermanAlignment {
    Cigar getCigar();
    int getAlignmentOffset();
}
