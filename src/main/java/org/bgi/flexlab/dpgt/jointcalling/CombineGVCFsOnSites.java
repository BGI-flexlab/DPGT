package org.bgi.flexlab.dpgt.jointcalling;

import java.lang.String;

public class CombineGVCFsOnSites {
    /**
     * JNI for combining gvcf of input region and variant sites
     * @param vcfpaths input vcf files
     * @param indexpaths input vcf index files
     * @param useLix 0: use vcf tbi/csi index, 1: use lix index
     * @param refpath  reference path
     * @param outpath  output file path of the combined gvcf
     * @param bytes variant site bitset as bytes
     * @param chrom chromosome
     * @param start 0-based start of the genome region
     * @param end   0-based end of the genome region
     */
    public native void Combine(String[] vcfpaths, String[] indexpaths, int useLix, String refpath,
        String outpath, byte[] bytes, String chrom, long start, long end);
}
