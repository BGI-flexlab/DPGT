package org.bgi.flexlab.gaea.tools.haplotypecaller.writer;

import java.io.File;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.util.ReadUtils;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.SAMFileHeader;

public final class SAMFileDestination extends HaplotypeBAMDestination {
    private final SAMFileGaeaReadWriter samWriter;

    /**
     * Create a new file destination.
     *
     * @param outFile file where output is written
     * @param createBamOutIndex true to create an index file for the bamout
     * @param createBamOutMD5 true to create an md5 file for the bamout
     * @param sourceHeader SAMFileHeader used to seed the output SAMFileHeader for this destination, must not be null
     * @param haplotypeReadGroupID  read group ID used when writing haplotypes as reads
     */
    public SAMFileDestination(
            final File outFile,
            final boolean createBamOutIndex,
            final boolean createBamOutMD5,
            final SAMFileHeader sourceHeader,
            final String haplotypeReadGroupID)
    {
        super(sourceHeader, haplotypeReadGroupID);
        samWriter = new SAMFileGaeaReadWriter(ReadUtils.createCommonSAMWriter(
                outFile,
                null,
                getBAMOutputHeader(), // use the header derived from the source header by HaplotypeBAMDestination
                false,
                createBamOutIndex,
                createBamOutMD5
        ));
    }

    /**
     * Close any resources associated with this destination.
     */
    @Override
    void close() { samWriter.close(); }

    /**
     * Write a read to the output file specified by this destination.
     *
     * @param read the read to write out, must not be null
     */
    @Override
    public void add(final GaeaSamRecord read) {
        Utils.nonNull(read, "read cannot be null");
        samWriter.addRead(read);
    }
}
