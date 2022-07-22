package org.bgi.flexlab.dpgt.jointcalling;

import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.bgi.flexlab.dpgt.utils.SimpleIntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.apache.spark.api.java.JavaFutureAction;
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
import java.io.File;
import java.io.InputStream;

import org.apache.log4j.PropertyConfigurator;
import org.apache.commons.io.FileUtils;

public class JointCallingSpark {
    private static final Logger logger = LoggerFactory.getLogger(JointCallingSpark.class);

    private static final String HEADER_DIR = "header";  // outdir/header
    private static final String GENOTYPE_HEADER = "genotype_header.vcf.gz";  // outdir/header/genotype_header.vcf.gz
    private static final String COMBINE_GVCFS_PREFIX = "combine."; // outdir/combine.idx
    private static final String GENOTYPE_GVCFS_PREFIX = "genotype."; // outdir/genotype.idx
    private static final String OUTPUT_NAME = "result.vcf.gz"; // outdir/result.vcf.gz

    private static final int GENOTYPE_PARTITION_COEFFICIENT = 4;

    static{
        try{
            System.loadLibrary("cdpgt");
        }catch(UnsatisfiedLinkError e){
            logger.error("Failed to load Native code library cdpgt.\n" + e);
            System.exit(1);
        }
    }

    public static void main(String[] args) throws Exception {
        // Initialize JNI C/C++ environment.
        NativeInitializer nativeInitializer = new NativeInitializer();
        nativeInitializer.apply();

        ClassLoader classloader = Thread.currentThread().getContextClassLoader();
        InputStream log4jPropertiestIs = classloader.getResourceAsStream("log4j.properties");
        PropertyConfigurator.configure(log4jPropertiestIs);
        JointCallingSparkOptions jcOptions = new JointCallingSparkOptions();
        jcOptions.parse(args);
        List<SimpleInterval> intervalsToTravers = jcOptions.getIntervalsToTravers();

        SparkConf conf = new SparkConf().setAppName("JointCallingSpark");
        if (jcOptions.uselocalMaster) conf.setMaster("local["+jcOptions.jobs+"]");
        JavaSparkContext sc = new JavaSparkContext(conf);

        JavaRDD<String> vcfpathsRDDPartitionByCombineParts = sc.textFile(jcOptions.input, jcOptions.numCombinePartitions);

        File headerDir = new File(jcOptions.output+"/"+HEADER_DIR);
        if (!headerDir.exists()) {
            headerDir.mkdirs();
        }
        combineVCFHeaders(vcfpathsRDDPartitionByCombineParts, headerDir.getAbsolutePath(), jcOptions.genotypeArguments);
        final String genotypeHeader = headerDir.getAbsolutePath() + "/" + GENOTYPE_HEADER;

        JavaRDD<String> vcfpathsRDDPartitionByJobs = sc.textFile(jcOptions.input, jcOptions.jobs);

        final String outputPath = jcOptions.output + "/" + OUTPUT_NAME;

        for (int i = 0; i < intervalsToTravers.size(); ++i) {
            SimpleInterval interval = intervalsToTravers.get(i);
            logger.info("Cycle {}/{}", i+1, intervalsToTravers.size());
            logger.info("Processing interval: {}", interval.toString());
            logger.info("Finding variant sites in {}", interval.toString());
            BitSet variantSiteSetData = vcfpathsRDDPartitionByJobs.
                mapPartitionsWithIndex(new VariantSiteFinderSparkFunc(interval.getContig(), interval.getStart()-1, interval.getEnd()-1), false).
                reduce((x, y) -> {x.or(y); return x;});

            if (variantSiteSetData.isEmpty()) {
                logger.info("skip interval: {}, because there is no variant site in it.", interval.toString());
                continue;
            }

            logger.info("Combining gvcfs in {}", interval.toString());
            File combineDir = new File(jcOptions.output+"/"+COMBINE_GVCFS_PREFIX+i);
            if (!combineDir.exists()) {
                combineDir.mkdirs();
            }
            final String combinePrefix = combineDir.getAbsolutePath()+"/"+COMBINE_GVCFS_PREFIX;

            ArrayList<SimpleInterval> subIntervals = SimpleIntervalUtils.splitIntervalByPartitions(interval,
                Math.max((int)Math.ceil(1.0*jcOptions.jobs/jcOptions.numCombinePartitions), 1));
            ArrayList<BitSet> subVariantBitSets = getSubBitSets(variantSiteSetData, interval, subIntervals);
            ArrayList<Broadcast<byte[]>> subVariantBitSetsBcs = new ArrayList<>();
            for (BitSet b: subVariantBitSets)
            {
                subVariantBitSetsBcs.add(sc.broadcast(b.toByteArray()));
            }

            ArrayList<JavaFutureAction<List<String>>> futureActions = new ArrayList<>();
            for (int j = 0; j < subIntervals.size(); ++j) {
                JavaFutureAction<List<String>> combinedGVCFsFuture = vcfpathsRDDPartitionByCombineParts.mapPartitionsWithIndex(new CombineGVCFsOnSitesSparkFunc(
                    jcOptions.reference, combinePrefix+j+"-", subIntervals.get(j), subVariantBitSetsBcs.get(j)), false)
                    .filter(x -> {return !x.equals("null");})
                    .collectAsync();  // collect async to run on multiple sub intervals at the same time
                futureActions.add(combinedGVCFsFuture);
            }
            ArrayList<List<String>> combinedGVCFsList = new ArrayList<>();
            for (int j = 0; j < futureActions.size(); ++j) {
                try {
                    combinedGVCFsList.add(futureActions.get(j).get());
                } catch (Exception e) {
                    logger.error("Failed to get combined gvcfs for interval: {}, {}", subIntervals.get(j), e.getMessage());
                    System.exit(1);
                }
            }

            logger.info("Genotyping gvcfs in {}", interval.toString());
            File genotypeDir = new File(jcOptions.output+"/"+GENOTYPE_GVCFS_PREFIX+i);
            if (!genotypeDir.exists()) {
                genotypeDir.mkdirs();
            }
            final String genotypePrefix = genotypeDir.getAbsolutePath()+"/"+GENOTYPE_GVCFS_PREFIX;
            ArrayList<String> genotypeGVCFsList = new ArrayList<>();
            for (int j = 0; j < subIntervals.size(); ++j) {
                ArrayList<SimpleInterval> windows = SimpleIntervalUtils.splitIntervalByPartitionsAndBitSet(
                    subIntervals.get(j), jcOptions.jobs * GENOTYPE_PARTITION_COEFFICIENT, subVariantBitSets.get(j));
                JavaRDD<SimpleInterval> windowsRDD = sc.parallelize(windows, windows.size());
                List<String> genotypeGVCFs = windowsRDD.
                    mapPartitionsWithIndex(new GVCFsSyncGenotyperSparkFunc(jcOptions.reference, combinedGVCFsList.get(j),
                        genotypeHeader, genotypePrefix+j+"-", jcOptions.dbsnp, jcOptions.genotypeArguments), false)
                    .filter(x -> {return !x.equals("null");})
                    .collect();
                genotypeGVCFsList.addAll(genotypeGVCFs);
            }
            
            logger.info("Concating genotyped vcfs in {} to result", interval.toString());
            JavaRDD<String> concateGVCFsRDD = sc.parallelize(genotypeGVCFsList, 1);
            concateGVCFsRDD.mapPartitionsWithIndex(new ConcatGenotypeGVCFsSparkFunc((i==0 ? genotypeHeader : null), outputPath), false).collect();

            if (jcOptions.deleteIntermediateResults) {
                FileUtils.deleteDirectory(combineDir);
                FileUtils.deleteDirectory(genotypeDir);
            }
        }

        if (jcOptions.deleteIntermediateResults) {
            FileUtils.deleteDirectory(headerDir);
        }

        sc.close();
    }

    private static VCFHeader combineVCFHeaders(JavaRDD<String> vcfpathsRDD, final String outdir, final GenotypeCalculationArgumentCollection genotypeArgs) {
        List<String> combinedHeaders = vcfpathsRDD.mapPartitionsWithIndex(new CombineVCFHeadersSparkFunc(outdir), false)
            .filter(x -> {return !x.equals("null");})
            .collect();

        // combine headers of each partition to generate combined header for all input vcfs
        VCFHeaderCombiner combiner = new VCFHeaderCombiner();
        String combinedHeader = outdir+"/header.vcf.gz";
        String[] combinedHeadersArray = new String[combinedHeaders.size()];
        combinedHeaders.toArray(combinedHeadersArray);
        combiner.Combine(combinedHeadersArray, combinedHeader);

        // make combined header for genotyping by add genotyping specific headers
        GVCFsSyncGenotyper genotyper = new GVCFsSyncGenotyper();
        VCFHeader combinedHeaderForGenotyping = genotyper.makeCombinedHeaderForGenotyping(combinedHeader, genotypeArgs);
        String combinedHeaderForGenotypingPath = outdir+"/"+GENOTYPE_HEADER;
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        builder.setOutputFile(combinedHeaderForGenotypingPath);
        builder.setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF);
        VariantContextWriter writer = builder.build();
        writer.writeHeader(combinedHeaderForGenotyping);
        writer.close();
        
        return combinedHeaderForGenotyping;
    }

    private static ArrayList<BitSet> getSubBitSets(final BitSet largeBitSet, final SimpleInterval largeInterval, final List<SimpleInterval> subIntervals)
    {
        ArrayList<BitSet> result = new ArrayList<>();
        for (SimpleInterval interval: subIntervals) {
            int i = interval.getStart() - largeInterval.getStart();
            int j = interval.getEnd() - largeInterval.getStart();
            result.add(largeBitSet.get(i, j));
        }
        return result;
    }
}
