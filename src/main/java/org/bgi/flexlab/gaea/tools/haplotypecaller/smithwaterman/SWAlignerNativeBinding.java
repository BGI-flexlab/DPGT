package org.bgi.flexlab.gaea.tools.haplotypecaller.smithwaterman;

import java.io.Closeable;
import java.io.File;

import org.bgi.flexlab.gaea.tools.haplotypecaller.smithwaterman.SmithWatermanJavaAligner.SWOverhangStrategy;

public interface SWAlignerNativeBinding extends  Closeable {

    /**
     * Perform a Smith-Waterman alignment and return the result
     *
     * @param ref reference sequence
     * @param alt alternate sequence
     * @return alignment result
     */
    SWNativeAlignerResult align(byte[] ref, byte[] alt, SWParameters parameters, SWOverhangStrategy overhangStrategy);


    /**
     * Subclasses may optionally implement close in order to release any native resources that they are holding.
     * Subclasses that rely on close to recover resources should fail with {@link IllegalStateException} if
     * {@link #align(byte[], byte[], SWParameters, SWOverhangStrategy)} is called after close.
     */
    @Override
    default void close() {
        //do nothing by default
    }
    
    boolean load(File tmpDir);
}
