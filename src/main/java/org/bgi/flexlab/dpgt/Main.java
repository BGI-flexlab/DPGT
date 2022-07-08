package org.bgi.flexlab.dpgt;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.slf4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.SparkConf;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.io.File;

import org.apache.commons.cli.*;


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

    public static void main(String[] args) throws Exception {
        Options options = new Options();
        options.addOption("i", "input", true, "Read a list of gvcf files from a file.");
        options.addOption("r", "reference", true, "Reference fasta file.");
        options.addOption("c", "combine-prefix", true, "Combine gvcfs file prefix.");
        options.addOption("g", "genotype-prefix", true, "Genotype gvcfs file prefix.");

        Parser parser = new PosixParser();
        HelpFormatter helpInfo = new HelpFormatter();
        CommandLine cmdLine = null;
        try {
            cmdLine = parser.parse(options, args);
        } catch (ParseException e) {
            helpInfo.printHelp("Options:", options);
            System.exit(1);
        }

        SparkConf conf = new SparkConf().setAppName("VariantSiteFinder").setMaster("local[10]");
        JavaSparkContext sc = new JavaSparkContext(conf);
        JavaRDD<String> vcfpathsRDD = sc.textFile(cmdLine.getOptionValue("input"), 10);
        
        File headerDir = new File("header");
        if (!headerDir.exists()) {
            headerDir.mkdirs();
        }
        combineVCFHeaders(vcfpathsRDD, "header");


        // SimpleInterval interval = new SimpleInterval("chr20", 1, 64444167);
        SimpleInterval interval = new SimpleInterval("chr20", 1, 5000000);
        logger.info("Finding variant sites...");
        BitSet variantSiteSetData = vcfpathsRDD.
            mapPartitionsWithIndex(new VariantSiteFinderSparkFunc(interval.getContig(), interval.getStart()-1, interval.getEnd()-1), false).
            reduce((x, y) -> {x.or(y); return x;});

        String refpath = cmdLine.getOptionValue("reference");
        String prefix = cmdLine.getOptionValue("combine-prefix");
        logger.info("Combining gvcfs on variant sites...");
        List<String> combinedGVCFs = vcfpathsRDD.mapPartitionsWithIndex(new CombineGVCFsOnSitesSparkFunc(
            refpath, prefix, variantSiteSetData.toByteArray(), interval.getContig(), interval.getStart()-1, interval.getEnd()-1), false).
            collect();

        ArrayList<SimpleInterval> windows = makeWindows(interval, 500000);
        JavaRDD<SimpleInterval> windowsRDD = sc.parallelize(windows, windows.size());
        List<String> genotypeGVCFs = windowsRDD.mapPartitionsWithIndex(new GVCFsSyncGenotyperSparkFunc(refpath, combinedGVCFs, cmdLine.getOptionValue("genotype-prefix")), false).
            collect();

        // for (int i = 0; i < windows.size(); ++i) {
        //     SimpleInterval testInterval = windows.get(i);
        //     logger.info("Genotyping interval: " + testInterval.toString());
        //     GVCFsSyncGenotyper genotyper = new GVCFsSyncGenotyper(refpath, combinedGVCFs, testInterval, args[3]+"."+ i +".vcf.gz");
        //     genotyper.run();
        //     genotyper.stop();
        // }
        
        logger.info("Done.");
        sc.close();

    }

    private static VCFHeader combineVCFHeaders(JavaRDD<String> vcfpathsRDD, final String outdir) {
        List<String> combinedHeaders = vcfpathsRDD.mapPartitionsWithIndex(new CombineVCFHeadersSparkFunc(outdir), false).collect();
        // combine headers of each partition to generate combined header for all input vcfs
        CombineVCFHeader combiner = new CombineVCFHeader();
        String combinedHeader = outdir+"/header.vcf.gz";
        combiner.call(combinedHeaders.iterator(), combinedHeader);

        GVCFsSyncGenotyper genotyper = new GVCFsSyncGenotyper();
        VCFHeader combinedHeaderForGenotyping = genotyper.makeCombinedHeaderForGenotyping(combinedHeader);
        String combinedHeaderForGenotypingPath = outdir+"/genotype_header.vcf.gz";
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        builder.setOutputFile(combinedHeaderForGenotypingPath);
        builder.setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF);
        builder.build().writeHeader(combinedHeaderForGenotyping);
        
        return combinedHeaderForGenotyping;
    }

    /**
     * make windows of specified size by splitting large interval
     * @param largeInterval large interval(window) to be splitted
     * @param size small interval(window) size
     * @return small windows
     */
    private static ArrayList<SimpleInterval> makeWindows(final SimpleInterval largeInterval, int size) {
        ArrayList<SimpleInterval> results = new ArrayList<>();
        if (largeInterval.size() < size) {
            results.add(largeInterval);
            return results;
        }
        int n = largeInterval.size() / size;
        for (int i = 0; i < n; ++i) {
            int s = largeInterval.getStart() + i * size;
            int e = s + size - 1;
            SimpleInterval window = new SimpleInterval(largeInterval.getContig(), s, e);
            results.add(window);
        }
        int lastEnd = results.get(n - 1).getEnd();
        if (lastEnd < largeInterval.getEnd()) {
            SimpleInterval window = new SimpleInterval(largeInterval.getContig(), lastEnd+1, largeInterval.getEnd());
            results.add(window);
        }
        return results;
    }

}
