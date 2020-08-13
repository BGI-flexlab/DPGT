package org.bgi.flexlab.gaea.tools.haplotypecaller.writer;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;

import htsjdk.samtools.SAMFileWriter;

public final class SAMFileGaeaReadWriter{

    private final SAMFileWriter samWriter;

    public SAMFileGaeaReadWriter( final SAMFileWriter samWriter ) {
        this.samWriter = samWriter;
    }

    public void addRead( GaeaSamRecord read ) {
        samWriter.addAlignment(read);
    }

    public void close() {
        samWriter.close();
    }
}
