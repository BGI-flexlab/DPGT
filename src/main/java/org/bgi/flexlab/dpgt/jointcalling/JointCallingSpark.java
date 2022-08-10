package org.bgi.flexlab.dpgt.jointcalling;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.SparkConf;
import java.util.ArrayList;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.List;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import org.bgi.flexlab.dpgt.utils.NativeLibraryLoader;

import org.apache.log4j.PropertyConfigurator;
import org.apache.commons.io.FileUtils;

public class JointCallingSpark {
    private static final Logger logger = LoggerFactory.getLogger(JointCallingSpark.class);

    private static final int PARTITION_COEFFICIENT = 4;

    static{
        NativeLibraryLoader.load();
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

        ConcatGenotypeGVCFsJob concatGVCFsJob0 = new ConcatGenotypeGVCFsJob(jcOptions);
        if (concatGVCFsJob0.isSuccess()) {
            logger.info("Joint Calling result {} is exists and previous job state was success. Please specify a new output directory or remove {}",
                jcOptions.getOutputVCFPath(), jcOptions.getOutputVCFPath());
            System.exit(0);
        }

        List<SimpleInterval> intervalsToTravers = jcOptions.getIntervalsToTravers();

        SparkConf conf = new SparkConf().setAppName("JointCallingSpark");
        if (jcOptions.uselocalMaster) conf.setMaster("local["+jcOptions.jobs+"]");
        if (conf.get("spark.master").startsWith("local")) jcOptions.uselocalMaster = true;
        NativeLibraryLoader.isLocal = jcOptions.uselocalMaster;
        JavaSparkContext sc = new JavaSparkContext(conf);

        JavaRDD<String> vcfpathsRDDPartitionByCombineParts = sc.textFile(addFilePrefixIfNeed(jcOptions.input), jcOptions.numCombinePartitions);

        CombineVCFHeaderJob combineVCFHeaderJob = new CombineVCFHeaderJob(vcfpathsRDDPartitionByCombineParts, jcOptions);
        final String genotypeHeader = combineVCFHeaderJob.run();

        JavaRDD<String> vcfpathsRDDPartitionByJobs = sc.textFile(addFilePrefixIfNeed(jcOptions.input), PARTITION_COEFFICIENT * jcOptions.jobs);

        ArrayList<File> allGenotypeDirs = new ArrayList<>();
        ArrayList<String> allGenotypeGVCFsList = new ArrayList<>();

        for (int i = 0; i < intervalsToTravers.size(); ++i) {
            SimpleInterval interval = intervalsToTravers.get(i);
            final String cycleStr = String.format("%d/%d", i+1, intervalsToTravers.size());
            logger.info("Cycle {}", cycleStr);
            logger.info("Processing interval: {} ({})", interval.toString(), cycleStr);
            // used for check if genotyping job is done for i-th iteration, if it is done we simply goto next iteration without
            // run or check variant finding job and combining gvcfs job
            GVCFsSyncGenotyperJob genotyperJob0 = new GVCFsSyncGenotyperJob(jcOptions, i);
            if (genotyperJob0.isSuccess()) {
                logger.info("Genotype gvcfs was success for interval: {} ({})", interval.toString(), cycleStr);
                allGenotypeGVCFsList.addAll(genotyperJob0.load());
                continue;
            }

            logger.info("Finding variant sites in {} ({})", interval.toString(), cycleStr);
            VariantFinderJob variantFinderJob = new VariantFinderJob(vcfpathsRDDPartitionByJobs, jcOptions, sc, i, interval);
            BitSet variantSiteSetData = variantFinderJob.run();
            
            if (variantSiteSetData == null || variantSiteSetData.isEmpty()) {
                logger.info("skip interval: {}, because there is no variant site in it ({})", interval.toString(), cycleStr);
                continue;
            }

            logger.info("Combining gvcfs in {} ({})", interval.toString(), cycleStr);
            CombineGVCFsOnSiteJob combineGVCFsOnSiteJob = new CombineGVCFsOnSiteJob(
                vcfpathsRDDPartitionByCombineParts, jcOptions, sc, i, interval, PARTITION_COEFFICIENT, variantSiteSetData);
            combineGVCFsOnSiteJob.run();

            logger.info("Genotyping gvcfs in {} ({})", interval.toString(), cycleStr);
            GVCFsSyncGenotyperJob genotyperJob = new GVCFsSyncGenotyperJob(jcOptions, sc, i,
                combineGVCFsOnSiteJob.subIntervals, combineGVCFsOnSiteJob.subVariantBitSets, combineGVCFsOnSiteJob.combinedGVCFsList,
                PARTITION_COEFFICIENT, genotypeHeader);
            List<String> genotypeGVCFsList = genotyperJob.run();
            allGenotypeGVCFsList.addAll(genotypeGVCFsList);
            allGenotypeDirs.add(genotyperJob.genotypeDir);

            if (jcOptions.deleteIntermediateResults) {
                FileUtils.deleteDirectory(variantFinderJob.variantSiteDir);
                FileUtils.deleteDirectory(combineGVCFsOnSiteJob.combineDir);
            }
        }

        logger.info("Concatenating genotyped vcfs to result");
        ConcatGenotypeGVCFsJob concatGVCFsJob = new ConcatGenotypeGVCFsJob(jcOptions, sc, allGenotypeGVCFsList, genotypeHeader, jcOptions.getOutputVCFPath());
        concatGVCFsJob.run();

        if (jcOptions.deleteIntermediateResults) {
            FileUtils.deleteDirectory(combineVCFHeaderJob.headerDir);
            for (final File d: allGenotypeDirs) {
                FileUtils.deleteDirectory(d);
            }
        }

        sc.close();
    }


    /**
     * add file:// prefix if need
     * @param input input file name
     * @return file name with file:// prefix
     */
    private static String addFilePrefixIfNeed(final String input) {
        if (!input.startsWith("file://")) {
            File inputFile = new File(input);
            try {
                return "file://" + inputFile.getCanonicalPath();
            } catch (IOException e) {
                logger.error(e.getMessage());
                System.exit(1);
            }
        }
        return input;
    }
}
