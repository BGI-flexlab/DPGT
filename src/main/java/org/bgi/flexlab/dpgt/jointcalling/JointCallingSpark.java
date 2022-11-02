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

        List<SimpleInterval> intervalsToTravers = jcOptions.getIntervalsToTravers();
        if (intervalsToTravers.size() < 1) {
            logger.error("Can not get intervals to travers!");
            System.exit(1);
        } 

        SparkConf conf = new SparkConf().setAppName("JointCallingSpark");
        if (jcOptions.uselocalMaster) conf.setMaster("local["+jcOptions.jobs+"]");
        if (conf.get("spark.master").startsWith("local")) jcOptions.uselocalMaster = true;
        NativeLibraryLoader.isLocal = jcOptions.uselocalMaster;
        JavaSparkContext sc = new JavaSparkContext(conf);

        JavaRDD<String> vcfpathsRDDPartitionByCombineParts = sc.textFile(addFilePrefixIfNeed(jcOptions.input), jcOptions.numCombinePartitions);

        CombineVCFHeaderJob combineVCFHeaderJob = new CombineVCFHeaderJob(vcfpathsRDDPartitionByCombineParts, jcOptions);
        final String genotypeHeader = combineVCFHeaderJob.run();

        JavaRDD<String> vcfpathsRDDPartitionByJobs = sc.textFile(addFilePrefixIfNeed(jcOptions.input), PARTITION_COEFFICIENT * jcOptions.jobs);

        int s = 0; // first interval that have not beed traversed
        int n = 0;
        for (int i = 0; i < intervalsToTravers.size(); ++i) {
            logger.info("Finding first interval to traverse ...");
            SimpleInterval interval = intervalsToTravers.get(i);
            final String cycleStr = String.format("%d/%d", i+1, intervalsToTravers.size());
            ConcatGenotypeGVCFsJob concatVcfJob0 = new ConcatGenotypeGVCFsJob(jcOptions, i);
            if (concatVcfJob0.isSuccess()) {
                ++n;
                logger.info("Genotype gvcfs was success for interval: {} ({})", interval.toString(), cycleStr);
            } else {
                s = i;
                logger.info("Start from interval: {} ({})", interval.toString(), cycleStr);
                break;
            }
        }

        if (n == intervalsToTravers.size()) {
            logger.info("Joint Calling result {} is exists and previous job state was success. Please specify a new output directory or remove {}",
                jcOptions.getOutputVCFPath(), jcOptions.getOutputVCFPath());
            System.exit(1);
        }

        ArrayList<File> allCombineDirs = new ArrayList<>();
        ArrayList<GVCFsSyncGenotyperJob> allGenotyperJobs = new ArrayList<>();
        ArrayList<File> allGenotypeDirs = new ArrayList<>();
        ConcatGenotypeGVCFsJob pre2ConcatVCFJob = null;

        for (int i = s; i < intervalsToTravers.size(); ++i) {
            SimpleInterval interval = intervalsToTravers.get(i);
            final String cycleStr = String.format("%d/%d", i+1, intervalsToTravers.size());
            logger.info("Cycle {}", cycleStr);
            logger.info("Processing interval: {} ({})", interval.toString(), cycleStr);
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
            allCombineDirs.add(combineGVCFsOnSiteJob.combineDir);

            if (jcOptions.deleteIntermediateResults) {
                FileUtils.deleteDirectory(variantFinderJob.variantSiteDir);
            }

            logger.info("Genotyping gvcfs in {} ({})", interval.toString(), cycleStr);
            GVCFsSyncGenotyperJob genotyperJob = new GVCFsSyncGenotyperJob(jcOptions, sc, i,
                combineGVCFsOnSiteJob.subIntervals, combineGVCFsOnSiteJob.subVariantBitSets, combineGVCFsOnSiteJob.combinedGVCFsList,
                PARTITION_COEFFICIENT, genotypeHeader);
            genotyperJob.submit();
            allGenotyperJobs.add(genotyperJob);
            allGenotypeDirs.add(genotyperJob.genotypeDir);

            if (i != s) {
                int j = i - s - 1;  // j >= 0
                GVCFsSyncGenotyperJob preGenotyperJob = allGenotyperJobs.get(j);
                List<String> preGenotypedVcfList = preGenotyperJob.get();  // -1 job
                ConcatGenotypeGVCFsJob preConcatVCFJob = null;
                if (i == 1) {
                    preConcatVCFJob = new ConcatGenotypeGVCFsJob(jcOptions, sc, preGenotypedVcfList, genotypeHeader, 0);
                } else {
                    preConcatVCFJob = new ConcatGenotypeGVCFsJob(jcOptions, sc, preGenotypedVcfList, null, i - 1);
                }
                if (pre2ConcatVCFJob != null) pre2ConcatVCFJob.get();  // -2 job
                preConcatVCFJob.submit();
                pre2ConcatVCFJob = preConcatVCFJob;

                // deleter intermediate results
                if (jcOptions.deleteIntermediateResults) {
                    FileUtils.deleteDirectory(allCombineDirs.get(j));  // delete -1 combine dir
                    if (j - 1 > -1) {
                        FileUtils.deleteDirectory(allGenotypeDirs.get(j-1));  // delete -2 genotype dir
                    }
                }
            }
        }

        {
            int j = intervalsToTravers.size() - s - 1;
            GVCFsSyncGenotyperJob preGenotyperJob = allGenotyperJobs.get(j);
            List<String> preGenotypedVcfList = preGenotyperJob.get();  // -1 job
            ConcatGenotypeGVCFsJob preConcatVCFJob = null;
            if (intervalsToTravers.size() == 1) {
                preConcatVCFJob = new ConcatGenotypeGVCFsJob(jcOptions, sc, preGenotypedVcfList, genotypeHeader, 0);
            } else {
                preConcatVCFJob = new ConcatGenotypeGVCFsJob(jcOptions, sc, preGenotypedVcfList, null, intervalsToTravers.size() - 1);
            }
            if (pre2ConcatVCFJob != null) pre2ConcatVCFJob.get();  // -2 job
            preConcatVCFJob.submit();
            preConcatVCFJob.get();

            // deleter intermediate results
            if (jcOptions.deleteIntermediateResults) {
                FileUtils.deleteDirectory(allCombineDirs.get(j));  // delete -1 combine dir
                if (j - 1 > -1) {
                    FileUtils.deleteDirectory(allGenotypeDirs.get(j-1));  // delete -2 genotype dir
                }
                FileUtils.deleteDirectory(allGenotypeDirs.get(j));   // delete -1 genotype dir
            }
        }

        if (jcOptions.deleteIntermediateResults) {
            FileUtils.deleteDirectory(combineVCFHeaderJob.headerDir);
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
