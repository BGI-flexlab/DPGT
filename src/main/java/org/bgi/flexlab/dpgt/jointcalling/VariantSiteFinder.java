package org.bgi.flexlab.dpgt.jointcalling;
import java.lang.String;

public class VariantSiteFinder {
    /**
     * JNI for finding variant site of input vcf files and input genome region
     * @param vcfpaths input vcf files
     * @param chrom chromosome
     * @param start 0-based start of the region
     * @param end   0-based end of the region
     * @return bitset as byte array
     */
    public native byte[] FindVariantSite(
        String[] vcfpaths, String chrom, long start, long end);
}
