package org.bgi.flexlab.dpgt.jointcalling;

import java.util.*;
import java.util.Iterator;
import org.apache.spark.api.java.function.Function2;
import org.bgi.flexlab.dpgt.utils.NativeLibraryLoader;


public class VariantSiteFinderSparkFunc implements Function2<Integer, Iterator<String>, Iterator<String>> {
    public String prefix;
    public String chrom;
    public int start;
    public int end;

    static {
        NativeLibraryLoader.load();
    }

    /**
     * variant site finder spark function
     * @param chrom chromosome
     * @param start 0-based start
     * @param end 0-based end
     */
    public VariantSiteFinderSparkFunc(final String prefix, final String chrom, int start, int end) {
        this.prefix = prefix;
        this.chrom = chrom;
        this.start = start;
        this.end = end;
    }
    @Override public Iterator<String> call(Integer idx, Iterator<String> vcfpathIter) {
        ArrayList<String> vcfpaths = new ArrayList<>();
        while(vcfpathIter.hasNext()) {
            vcfpaths.add(vcfpathIter.next());
        }
        if (vcfpaths.isEmpty()) {
            // no vcf in this part, return a null output path
            ArrayList<String> result = new ArrayList<>();
            result.add("null");
            return result.iterator();
        }
        String[] vcfpathsArray  = new String[vcfpaths.size()];
        vcfpaths.toArray(vcfpathsArray);
        String outpath = prefix + idx + JointCallingSparkConsts.VARIANT_SITE_SUFFIX;
        VariantSiteFinder vf = new VariantSiteFinder();
        // call native c++ function to find variant site, return variant site bitset as byte array
        vf.FindVariantSite(vcfpathsArray, outpath, this.chrom, (long)this.start, (long)this.end);
        ArrayList<String> result = new ArrayList<>();
        result.add(outpath);
        return result.iterator();
    }
}
