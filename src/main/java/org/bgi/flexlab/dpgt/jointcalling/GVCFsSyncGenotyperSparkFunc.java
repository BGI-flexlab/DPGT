package org.bgi.flexlab.dpgt.jointcalling;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
import org.apache.spark.api.java.function.Function2;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;


public class GVCFsSyncGenotyperSparkFunc implements Function2<Integer, Iterator<SimpleInterval>, Iterator<String>> {
    public String refpath;
    public List<String> vcfpaths;
    public String vcfHeaderPath;
    public String prefix;
    public String dbsnpPath;
    public GenotypeCalculationArgumentCollection genotypeArgs;
    
    public GVCFsSyncGenotyperSparkFunc(final String refpath, final List<String> vcfpaths, final String vcfHeaderPath,
        final String prefix, final String dbsnpPath, final GenotypeCalculationArgumentCollection genotypeArgs) {
        this.refpath = refpath;
        this.vcfpaths = vcfpaths;
        this.vcfHeaderPath = vcfHeaderPath;
        this.prefix = prefix;
        this.dbsnpPath = dbsnpPath;
        this.genotypeArgs = genotypeArgs;
    }

    @Override public Iterator<String> call(Integer idx, Iterator<SimpleInterval> intervalIter) {
        if (intervalIter.hasNext()) {
            SimpleInterval interval = intervalIter.next();
            String outpath = prefix + idx + ".vcf.gz";
            GVCFsSyncGenotyper genotyper = new GVCFsSyncGenotyper(refpath, vcfpaths, vcfHeaderPath, interval, outpath, dbsnpPath, genotypeArgs);
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
