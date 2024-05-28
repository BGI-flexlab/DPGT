package org.bgi.flexlab.dpgt.jointcalling;

import java.util.*;

import org.apache.spark.api.java.function.Function2;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.bgi.flexlab.dpgt.utils.NativeLibraryLoader;


public class CombineGVCFsOnSitesSparkFunc implements Function2<Integer, Iterator<String>, Iterator<String>> {
    public String refpath;
    public String prefix;
    public int useLix = 0;
    public SimpleInterval interval;
    Broadcast<byte[]> variantSiteSetBytesBc;

    static {
        NativeLibraryLoader.load();
    }

    public CombineGVCFsOnSitesSparkFunc(final String refpath, final String prefix,
        int useLix, final SimpleInterval interval, Broadcast<byte[]> variantSiteSetBytesBc)
    {
        this.refpath = refpath;
        this.prefix = prefix;
        this.useLix = useLix;
        this.interval = interval;
        this.variantSiteSetBytesBc = variantSiteSetBytesBc;
    }

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
            ArrayList<String> returnValue=new ArrayList<>();
            returnValue.add("null");
            return returnValue.iterator();
        }
        String[] vcfpathsArray  = new String[vcfpaths.size()];
        String[] vcfindicesArray = new String[vcfindices.size()];
        vcfpaths.toArray(vcfpathsArray);
        vcfindices.toArray(vcfindicesArray);
        CombineGVCFsOnSites combiner = new CombineGVCFsOnSites();
        String outpath = prefix + idx + ".vcf.gz";
        // call native c++ function to combine gvcfs
        combiner.Combine(vcfpathsArray, vcfindicesArray, useLix, refpath, outpath, variantSiteSetBytesBc.value(), interval.getContig(), interval.getStart()-1, interval.getEnd()-1);
        ArrayList<String> result=new ArrayList<>();
        result.add(outpath);
        return result.iterator();
    }
}
