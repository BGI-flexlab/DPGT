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
        logger.info("Finding first interval to traverse ...");
        for (int i = 0; i < intervalsToTravers.size(); ++i) {
            SimpleInterval interval = intervalsToTravers.get(i);
            final String cycleStr = String.format("%d/%d", i+1, intervalsToTravers.size());
            ConcatGenotypeGVCFsJob concatVcfJob0 = new ConcatGenotypeGVCFsJob(jcOptions, i);
            if (concatVcfJob0.isSuccess()) {
                ++n;
                logger.info("Genotyping gvcfs was success and results were written to output file for interval: {} ({})",
                    interval.toString(), cycleStr);
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
        ConcatGenotypeGVCFsJob pre2ConcatVCFJob = new ConcatGenotypeGVCFsJob(); // construct a null job

        for (int i = s; i < intervalsToTravers.size(); ++i) {
            SimpleInterval interval = intervalsToTravers.get(i);
            final String cycleStr = String.format("%d/%d", i+1, intervalsToTravers.size());
            logger.info("Cycle {}", cycleStr);
            logger.info("Processing interval: {} ({})", interval.toString(), cycleStr);

            GVCFsSyncGenotyperJob genotyperJob0 = new GVCFsSyncGenotyperJob(jcOptions, i);
            if (genotyperJob0.isSuccess()) {
                allCombineDirs.add(null);
                allGenotyperJobs.add(genotyperJob0);
                allGenotypeDirs.add(genotyperJob0.genotypeDir);
                logger.info("skip, because genotyping gvcfs was success for interval: {} ({})",
                    interval.toString(), cycleStr);
                handleJobDependency(i, s, allGenotyperJobs, jcOptions, sc, genotypeHeader, pre2ConcatVCFJob, allCombineDirs, allGenotypeDirs);
                continue;
            }

            logger.info("Finding variant sites in {} ({})", interval.toString(), cycleStr);
            VariantFinderJob variantFinderJob = new VariantFinderJob(vcfpathsRDDPartitionByJobs, jcOptions, sc, i, interval);
            BitSet variantSiteSetData = variantFinderJob.run();
            
            if (variantSiteSetData == null || variantSiteSetData.isEmpty()) {
                allCombineDirs.add(null);
                allGenotyperJobs.add(null);
                allGenotypeDirs.add(null);
                logger.info("skip interval: {}, because there is no variant site in it ({})", interval.toString(), cycleStr);
                handleJobDependency(i, s, allGenotyperJobs, jcOptions, sc, genotypeHeader, pre2ConcatVCFJob, allCombineDirs, allGenotypeDirs);
                continue;
            }

            logger.info("Combining gvcfs in {} ({})", interval.toString(), cycleStr);
            CombineGVCFsOnSiteJob combineGVCFsOnSiteJob = new CombineGVCFsOnSiteJob(
                vcfpathsRDDPartitionByCombineParts, jcOptions, sc, i, interval, PARTITION_COEFFICIENT, variantSiteSetData);
            combineGVCFsOnSiteJob.run();
            allCombineDirs.add(combineGVCFsOnSiteJob.combineDir);

            if (jcOptions.deleteIntermediateResults) {
                safeDeleteDirectory(variantFinderJob.variantSiteDir);
            }

            logger.info("Genotyping gvcfs in {} ({})", interval.toString(), cycleStr);
            GVCFsSyncGenotyperJob genotyperJob = new GVCFsSyncGenotyperJob(jcOptions, sc, i,
                combineGVCFsOnSiteJob.subIntervals, combineGVCFsOnSiteJob.subVariantBitSets, combineGVCFsOnSiteJob.combinedGVCFsList,
                PARTITION_COEFFICIENT, genotypeHeader);
            genotyperJob.submit();
            allGenotyperJobs.add(genotyperJob);
            allGenotypeDirs.add(genotyperJob.genotypeDir);

            pre2ConcatVCFJob = handleJobDependency(i, s, allGenotyperJobs, jcOptions, sc, genotypeHeader, pre2ConcatVCFJob, allCombineDirs, allGenotypeDirs);
        }

        pre2ConcatVCFJob = handleJobDependency(intervalsToTravers.size(), s, allGenotyperJobs, jcOptions, sc, genotypeHeader, pre2ConcatVCFJob, allCombineDirs, allGenotypeDirs);

        if (!pre2ConcatVCFJob.isNull()) pre2ConcatVCFJob.get();
        int j = intervalsToTravers.size() - s - 1;
        safeDeleteDirectory(allGenotypeDirs.get(j));

        if (jcOptions.deleteIntermediateResults) {
            safeDeleteDirectory(combineVCFHeaderJob.headerDir);
        }

        sc.close();
    }


    private static ConcatGenotypeGVCFsJob handleJobDependency(int i, int s, 
        ArrayList<GVCFsSyncGenotyperJob> allGenotyperJobs,
        JointCallingSparkOptions jcOptions,
        JavaSparkContext sc,
        String genotypeHeader,
        ConcatGenotypeGVCFsJob pre2ConcatVCFJob,
        List<File> allCombineDirs,
        List<File> allGenotypeDirs
        ) throws IOException
    {
        ConcatGenotypeGVCFsJob preConcatVCFJob = null;
        if (i != s) {
            int j = i - s - 1;  // j >= 0
            GVCFsSyncGenotyperJob preGenotyperJob = allGenotyperJobs.get(j);
            if (preGenotyperJob != null) {
                List<String> preGenotypedVcfList = preGenotyperJob.get();  // -1 job
                if (i == 1) {
                    preConcatVCFJob = new ConcatGenotypeGVCFsJob(jcOptions, sc, preGenotypedVcfList, genotypeHeader, 0);
                } else {
                    preConcatVCFJob = new ConcatGenotypeGVCFsJob(jcOptions, sc, preGenotypedVcfList, null, i - 1);
                }
                if (!pre2ConcatVCFJob.isNull()) pre2ConcatVCFJob.get();  // -2 job
                preConcatVCFJob.submit();
                pre2ConcatVCFJob = preConcatVCFJob;
            } else {
                if (i == 1) {
                    ArrayList<String> preGenotypedVcfList = new ArrayList<>();
                    preGenotypedVcfList.add(null);
                    preConcatVCFJob = new ConcatGenotypeGVCFsJob(jcOptions, sc, preGenotypedVcfList, genotypeHeader, 0);
                }
                if (preConcatVCFJob != null) {
                    pre2ConcatVCFJob = preConcatVCFJob;
                } else {
                    pre2ConcatVCFJob.setToNull();
                }
            }
            // deleter intermediate results
            if (jcOptions.deleteIntermediateResults) {
                safeDeleteDirectory(allCombineDirs.get(j));  // delete -1 combine dir
                if (j - 1 > -1) {
                    safeDeleteDirectory(allGenotypeDirs.get(j-1));  // delete -2 genotype dir
                }
            }
        }
        return pre2ConcatVCFJob;
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

    private static void safeDeleteDirectory(File dir) throws IOException {
        if (dir != null && dir.exists()) {
            FileUtils.deleteDirectory(dir);
        }
    }
}
