package org.bgi.flexlab.dpgt;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.slf4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.SparkConf;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.List;


public class Main {
    private static final Logger logger = LoggerFactory.getLogger(Main.class);

    static{
        try{
            System.loadLibrary("cdpgt");
        }catch(UnsatisfiedLinkError e){
            logger.error("Native code library failed to load.\n" + e);
            System.exit(1);
        }
    }

    public static void main(String[] args) {
        SparkConf conf = new SparkConf().setAppName("VariantSiteFinder").setMaster("local[10]");
        JavaSparkContext sc = new JavaSparkContext(conf);
        JavaRDD<String> vcfpathsRDD = sc.textFile(args[0], 10);
        String chrom = "chr20";
        long start = 0;
        long end = 500000;
        long size = end - start + 1;
        logger.info("Finding variant sites...");
        BitSet variantSiteSetData = vcfpathsRDD.
            mapPartitionsWithIndex(new VariantSiteFinderSparkFunc(chrom, start, end), false).
            reduce((x, y) -> {x.or(y); return x;});

        String refpath = args[1];
        String prefix = args[2];
        logger.info("Combining gvcfs on variant sites...");
        List<String> combinedGVCFs = vcfpathsRDD.mapPartitionsWithIndex(new CombineGVCFsOnSitesSparkFunc(
            refpath, prefix, variantSiteSetData.toByteArray(), chrom, start, end), false).
            collect();
        SimpleInterval interval = new SimpleInterval(chrom, (int)start+1, (int)end+1);
        GVCFsSyncGenotyper genotyper = new GVCFsSyncGenotyper(refpath, combinedGVCFs, interval);
        
        sc.close();
        logger.info("Done.");
    }
}
