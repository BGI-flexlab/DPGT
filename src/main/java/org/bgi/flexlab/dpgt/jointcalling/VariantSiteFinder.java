package org.bgi.flexlab.dpgt.jointcalling;
import java.lang.String;

public class VariantSiteFinder {
    /**
     * JNI for finding variant site of input vcf files and input genome region
     * @param vcfpaths input vcf files
     * @param indexpaths input vcf index files
     * @param useLix 0: use vcf tbi/csi index, 1: use lix index
     * @param outpath output variant site bit set bytes in binary format
     * @param chrom chromosome
     * @param start 0-based start of the region
     * @param end   0-based end of the region
     * @return
     */
    public native void FindVariantSite(
        String[] vcfpaths, String[] indexpaths, int useLix,
        String outpath, String chrom, long start, long end);
}
