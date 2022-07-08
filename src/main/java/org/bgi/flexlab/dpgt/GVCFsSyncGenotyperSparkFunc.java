package org.bgi.flexlab.dpgt;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
import org.apache.spark.api.java.function.Function2;
import org.broadinstitute.hellbender.utils.SimpleInterval;


public class GVCFsSyncGenotyperSparkFunc implements Function2<Integer, Iterator<SimpleInterval>, Iterator<String>> {
    public String refpath;
    public List<String> vcfpaths;
    public String vcfHeaderPath;
    public String outdir;
    
    public GVCFsSyncGenotyperSparkFunc(final String refpath, final List<String> vcfpaths, final String vcfHeaderPath, final String outdir) {
        this.refpath = refpath;
        this.vcfpaths = vcfpaths;
        this.vcfHeaderPath = vcfHeaderPath;
        this.outdir = outdir;
    }

    @Override public Iterator<String> call(Integer idx, Iterator<SimpleInterval> intervalIter) {
        if (intervalIter.hasNext()) {
            SimpleInterval interval = intervalIter.next();
            String outpath = outdir + "/genotypeGVCFs." + idx + ".vcf.gz";
            GVCFsSyncGenotyper genotyper = new GVCFsSyncGenotyper(refpath, vcfpaths, interval, outpath);
            genotyper.run();
            genotyper.stop();
            ArrayList<String> result=new ArrayList<>();
            result.add(outpath);
            return result.iterator();
        }
        ArrayList<String> result=new ArrayList<>();
        result.add("null");
        return result.iterator();
    }
}
