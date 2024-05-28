package org.bgi.flexlab.dpgt.jointcalling;

import java.util.*;
import java.util.Iterator;
import org.apache.spark.api.java.function.Function2;
import org.bgi.flexlab.dpgt.utils.NativeLibraryLoader;


public class VariantSiteFinderSparkFunc implements Function2<Integer, Iterator<String>, Iterator<String>> {
    public String prefix;
    public int useLix = 0;
    public String chrom;
    public int start;
    public int end;

    static {
        NativeLibraryLoader.load();
    }

    /**
     * variant site finder spark function
     * @param prefix output prefix
     * @param use_lix 0 or 1, 1 means use lix index, 0 means use vcf tbi/csi index
     * @param chrom chromosome
     * @param start 0-based start
     * @param end 0-based end
     */
    public VariantSiteFinderSparkFunc(final String prefix, int useLix, final String chrom, int start, int end) {
        this.prefix = prefix;
        this.useLix = useLix;
        this.chrom = chrom;
        this.start = start;
        this.end = end;
    }

    /**
     * call native c++ function to find variant site
     * @param idx partition index
     * @param vcfPairIter iterator of vcf file and vcf index pairs(commas separated vcf file and vcf index) 
     */
    @Override public Iterator<String> call(Integer idx, Iterator<String> vcfPairIter) {
        ArrayList<String> vcfpaths = new ArrayList<>();
        ArrayList<String> vcfindices = new ArrayList<>();
        while(vcfPairIter.hasNext()) {
            String pair = vcfPairIter.next();
            if (pair == null || pair.isEmpty()) {
                continue;
            }
            if (pair.contains(",")) {
                String[] pairArray = pair.split(",");
                vcfpaths.add(pairArray[0]);
                vcfindices.add(pairArray[1]);
            } else {
                vcfpaths.add(pair);
            }
        }
        if (vcfpaths.isEmpty()) {
            // no vcf in this part, return a null output path
            ArrayList<String> result = new ArrayList<>();
            result.add("null");
            return result.iterator();
        }
        String[] vcfpathsArray  = new String[vcfpaths.size()];
        String[] vcfindicesArray = new String[vcfindices.size()];
        vcfpaths.toArray(vcfpathsArray);
        vcfindices.toArray(vcfindicesArray);
        String outpath = prefix + idx + JointCallingSparkConsts.VARIANT_SITE_SUFFIX;
        VariantSiteFinder vf = new VariantSiteFinder();
        // call native c++ function to find variant site, return variant site bitset as byte array
        vf.FindVariantSite(vcfpathsArray, vcfindicesArray, this.useLix, outpath, this.chrom, (long)this.start, (long)this.end);
        ArrayList<String> result = new ArrayList<>();
        result.add(outpath);
        return result.iterator();
    }
}
