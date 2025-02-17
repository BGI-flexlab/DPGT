/**
This file is part of DPGT.
Copyright (C) 2022 BGI.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// License End
package org.bgi.flexlab.dpgt.jointcalling;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.SparkConf;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.List;
import java.util.TreeMap;
import java.util.concurrent.TimeUnit;
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

        JavaRDD<String> vcfpathsRDDPartitionByJobs = sc.textFile(addFilePrefixIfNeed(jcOptions.input), jcOptions.jobs);

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

        TreeMap<Integer, GVCFsSyncGenotyperJob> allGenotyperJobs = new TreeMap<>();
        TreeMap<Integer, ConcatGenotypeGVCFsJob> allConcatVcfJobs = new TreeMap<>();

        for (int i = s; i < intervalsToTravers.size(); ++i) {
            SimpleInterval interval = intervalsToTravers.get(i);
            final String cycleStr = String.format("%d/%d", i+1, intervalsToTravers.size());
            logger.info("Cycle {}", cycleStr);
            logger.info("Processing interval: {} ({})", interval.toString(), cycleStr);

            ConcatGenotypeGVCFsJob concatVcfJob0 = new ConcatGenotypeGVCFsJob(jcOptions, i);
            if (concatVcfJob0.isSuccess()) {
                logger.info("skip, because genotyping gvcfs results were written to output file for interval: {} ({})",
                    interval.toString(), cycleStr);
                continue;
            }

            GVCFsSyncGenotyperJob genotyperJob0 = new GVCFsSyncGenotyperJob(jcOptions, i);
            if (genotyperJob0.isSuccess()) {
                allGenotyperJobs.put(i, genotyperJob0);
                logger.info("skip, because genotyping gvcfs was success for interval: {} ({})",
                    interval.toString(), cycleStr);
                continue;
            }

            logger.info("Finding variant sites in {} ({})", interval.toString(), cycleStr);
            VariantFinderJob variantFinderJob = new VariantFinderJob(vcfpathsRDDPartitionByJobs, jcOptions, sc, i, interval);
            BitSet variantSiteSetData = variantFinderJob.run();
            updateGenotypeAndConcatVcfJobs(allGenotyperJobs, allConcatVcfJobs, intervalsToTravers, jcOptions, genotypeHeader, sc, JointCallingSparkConsts.GET_FUTURE_TIMEOUT);
            
            if (variantSiteSetData == null || variantSiteSetData.isEmpty()) {
                logger.info("skip interval: {}, because there is no variant site in it ({})", interval.toString(), cycleStr);
                continue;
            }

            logger.info("Combining gvcfs in {} ({})", interval.toString(), cycleStr);
            CombineGVCFsOnSiteJob combineGVCFsOnSiteJob = new CombineGVCFsOnSiteJob(
                vcfpathsRDDPartitionByCombineParts, jcOptions, sc, i, interval, PARTITION_COEFFICIENT, variantSiteSetData);
            combineGVCFsOnSiteJob.run();
            updateGenotypeAndConcatVcfJobs(allGenotyperJobs, allConcatVcfJobs, intervalsToTravers, jcOptions, genotypeHeader, sc, JointCallingSparkConsts.GET_FUTURE_TIMEOUT);

            logger.info("Genotyping gvcfs in {} ({})", interval.toString(), cycleStr);
            GVCFsSyncGenotyperJob genotyperJob = new GVCFsSyncGenotyperJob(jcOptions, sc, i,
                combineGVCFsOnSiteJob.subIntervals, combineGVCFsOnSiteJob.subVariantBitSets, combineGVCFsOnSiteJob.combinedGVCFsList,
                PARTITION_COEFFICIENT, genotypeHeader);
            genotyperJob.addToRemoveList(combineGVCFsOnSiteJob.combineDir);
            genotyperJob.submit();
            allGenotyperJobs.put(i, genotyperJob);

            updateGenotypeAndConcatVcfJobs(allGenotyperJobs, allConcatVcfJobs, intervalsToTravers, jcOptions, genotypeHeader, sc, JointCallingSparkConsts.GET_FUTURE_TIMEOUT);
        }

        updateGenotypeAndConcatVcfJobs(allGenotyperJobs, allConcatVcfJobs, intervalsToTravers, jcOptions, genotypeHeader, sc, 0);

        for (Integer k: allConcatVcfJobs.keySet()) {
            ConcatGenotypeGVCFsJob concatGenotypeGVCFsJob = allConcatVcfJobs.get(k);
            concatGenotypeGVCFsJob.get();
        }

        if (jcOptions.deleteIntermediateResults) {
            safeDeleteDirectory(combineVCFHeaderJob.headerDir);
            for (int i = 0; i < intervalsToTravers.size(); ++i) {
                File variantSiteDir = new File(jcOptions.output+"/"+JointCallingSparkConsts.VARIANT_SITE_PREFIX+i);
                safeDeleteDirectory(variantSiteDir);
            }
        }

        sc.close();
    }


    private static void updateGenotypeAndConcatVcfJobs(TreeMap<Integer, GVCFsSyncGenotyperJob> allGenotyperJobs,
        TreeMap<Integer, ConcatGenotypeGVCFsJob> allConcatVcfJobs, List<SimpleInterval> intervalsToTravers,
        JointCallingSparkOptions jcOptions, String genotypeHeader, JavaSparkContext sc, long timeout)
    {
        for (Integer k: allGenotyperJobs.keySet()) {
            GVCFsSyncGenotyperJob gtJob = allGenotyperJobs.get(k);
            if (gtJob != null) {
                List<String> gtList = null;
                if (timeout > 0) {
                    gtList = gtJob.get(timeout, TimeUnit.MILLISECONDS);
                } else {
                    gtList = gtJob.get();
                }
                if (gtList != null) {
                    if (allConcatVcfJobs.get(k) == null) {
                        ConcatGenotypeGVCFsJob concatGenotypeGVCFsJob = new ConcatGenotypeGVCFsJob(
                            jcOptions, sc, gtList, genotypeHeader, k, intervalsToTravers.get(k));
                        concatGenotypeGVCFsJob.addToRemoveList(gtJob.genotypeDir);
                        concatGenotypeGVCFsJob.submit();
                        allConcatVcfJobs.put(k, concatGenotypeGVCFsJob);
                    } else {
                        if (timeout > 0) {
                            allConcatVcfJobs.get(k).get(timeout, TimeUnit.MILLISECONDS);
                        } else {
                            allConcatVcfJobs.get(k).get();
                        }
                    }
                }
            }
        }
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

    private static void safeDeleteDirectory(File dir) {
        if (dir != null && dir.exists()) {
            try {
                FileUtils.deleteDirectory(dir);
            } catch (Exception e) {
                // in some case, we may failed to delete dir because it is using by other process.
                logger.warn("Unable to delete directory {}, {}.", dir.getAbsolutePath(), e.getMessage());
            }
        }
    }
}
