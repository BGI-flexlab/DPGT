package org.bgi.flexlab.dpgt.jointcalling;

import java.util.*;
import java.util.Iterator;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.utils.SimpleInterval;


public class CombineGVCFsOnSitesSparkFunc implements Function2<Integer, Iterator<String>, Iterator<String>> {
    public String refpath;
    public String prefix;
    public SimpleInterval interval;
    Broadcast<byte[]> variantSiteSetBytesBc;
    public CombineGVCFsOnSitesSparkFunc(final String refpath,
        final String prefix, final SimpleInterval interval, Broadcast<byte[]> variantSiteSetBytesBc)
    {
        this.refpath = refpath;
        this.prefix = prefix;
        this.interval = interval;
        this.variantSiteSetBytesBc = variantSiteSetBytesBc;
    }

    @Override public Iterator<String> call(Integer idx, Iterator<String> vcfpathIter) {
        ArrayList<String> vcfpaths = new ArrayList<>();
        while(vcfpathIter.hasNext()) {
            vcfpaths.add(vcfpathIter.next());
        }
        String[] vcfpathsArray  = new String[vcfpaths.size()];
        vcfpaths.toArray(vcfpathsArray);
        CombineGVCFsOnSites combiner = new CombineGVCFsOnSites();
        String outpath = prefix + idx + ".vcf.gz";
        // call native c++ function to combine gvcfs
        combiner.Combine(vcfpathsArray, refpath, outpath, variantSiteSetBytesBc.value(), interval.getContig(), interval.getStart()-1, interval.getEnd()-1);
        ArrayList<String> result=new ArrayList<>();
        result.add(outpath);
        return result.iterator();
    }
}
