package org.bgi.flexlab.dpgt.jointcalling;

import java.util.*;
import java.util.Iterator;
import java.util.BitSet;
import org.apache.spark.api.java.function.Function2;


public class VariantSiteFinderSparkFunc implements Function2<Integer, Iterator<String>, Iterator<BitSet>> {
    public String chrom;
    public int start;
    public int end;
    /**
     * variant site finder spark function
     * @param chrom chromosome
     * @param start 0-based start
     * @param end 0-based end
     */
    public VariantSiteFinderSparkFunc(String chrom, int start, int end) {
        this.chrom = chrom;
        this.start = start;
        this.end = end;
    }
    @Override public Iterator<BitSet> call(Integer idx, Iterator<String> vcfpathIter) {
        ArrayList<String> vcfpaths = new ArrayList<>();
        while(vcfpathIter.hasNext()) {
            vcfpaths.add(vcfpathIter.next());
        }
        String[] vcfpathsArray  = new String[vcfpaths.size()];
        vcfpaths.toArray(vcfpathsArray);
        VariantSiteFinder vf = new VariantSiteFinder();
        // call native c++ function to find variant site, return variant site bitset as byte array
        byte[] res = vf.FindVariantSite(vcfpathsArray, this.chrom, (long)this.start, (long)this.end);
        BitSet variantSiteSet = BitSet.valueOf(res);
        ArrayList<BitSet> result = new ArrayList<>();
        result.add(variantSiteSet);
        return result.iterator();
    }
}
