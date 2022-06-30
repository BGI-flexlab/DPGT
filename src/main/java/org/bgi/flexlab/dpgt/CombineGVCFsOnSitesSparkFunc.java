package org.bgi.flexlab.dpgt;

import java.util.*;
import java.util.Iterator;
import org.apache.spark.api.java.function.Function2;


public class CombineGVCFsOnSitesSparkFunc implements Function2<Integer, Iterator<String>, Iterator<String>> {
    public String refpath;
    public String prefix;
    public byte[] bytes;
    public String chrom;
    public long start;
    public long end;
    public CombineGVCFsOnSitesSparkFunc(String refpath,
        String prefix, byte[] bytes, String chrom, long start, long end)
    {
        this.refpath = refpath;
        this.prefix = prefix;
        this.bytes = bytes;
        this.chrom = chrom;
        this.start = start;
        this.end = end;
    }

    @Override public Iterator<String> call(Integer idx, Iterator<String> vcfpathIter) {
        ArrayList<String> vcfpaths = new ArrayList<>();
        while(vcfpathIter.hasNext()) {
            vcfpaths.add(vcfpathIter.next());
        }
        String[] vcfpathsArray  = new String[vcfpaths.size()];
        vcfpaths.toArray(vcfpathsArray);
        CombineGVCFsOnSites combiner = new CombineGVCFsOnSites();
        String outpath = prefix + "." + idx + ".vcf.gz";
        // call native c++ function to combine gvcfs
        combiner.Combine(vcfpathsArray, refpath, outpath, bytes, chrom, start, end);
        ArrayList<String> result=new ArrayList<>();
        result.add(outpath);  // make a nonsense return value
        return result.iterator();
    }
}
