package org.bgi.flexlab.dpgt.jointcalling;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.bgi.flexlab.dpgt.utils.DPGTJobAsync;
import org.bgi.flexlab.dpgt.utils.DPGTJobState;
import org.bgi.flexlab.dpgt.utils.SimpleIntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.JavaFutureAction;

public class GVCFsSyncGenotyperJob extends DPGTJobAsync<List<String>, List<String>> {
    private static final Logger logger = LoggerFactory.getLogger(GVCFsSyncGenotyperJob.class);
    private JointCallingSparkOptions jcOptions;
    private JavaSparkContext sc;
    private int idx;
    private List<SimpleInterval> subIntervals;
    private List<BitSet> subVariantBitSets;
    private List<List<String>> combinedGVCFsList;
    private int partitionCoeff;
    private String genotypeHeader;
    public File genotypeDir;
    public GVCFsSyncGenotyperJob(final JointCallingSparkOptions jcOptions, final JavaSparkContext sc, final int idx,
        final List<SimpleInterval> subIntervals, final List<BitSet> subVariantBitSets, final List<List<String>> combinedGVCFsList,
        final int partitionCoeff, final String genotypeHeader)
    {
        this.jcOptions = jcOptions;
        this.sc = sc;
        this.idx = idx;
        this.subIntervals = subIntervals;
        this.subVariantBitSets = subVariantBitSets;
        this.combinedGVCFsList = combinedGVCFsList;
        this.partitionCoeff = partitionCoeff;
        this.genotypeHeader = genotypeHeader;
        this.genotypeDir = new File(jcOptions.output+"/"+JointCallingSparkConsts.GENOTYPE_GVCFS_PREFIX+this.idx);
        if (!this.genotypeDir.exists()) {
            this.genotypeDir.mkdirs();
        }
        this.stateFile = this.jcOptions.output + "/" + JointCallingSparkConsts.JOB_STATE + "/" + JointCallingSparkConsts.GENOTYPE_GVCFS_PREFIX + this.idx + ".json";
    }

    /**
     * this constructor is used for construct a job for check job state(isSuccess), DO NOT run it
     * @param jcOptions
     * @param idx
     */
    public GVCFsSyncGenotyperJob(final JointCallingSparkOptions jcOptions, final int idx) {
        this.jcOptions = jcOptions;
        this.idx = idx;
        this.genotypeDir = new File(jcOptions.output+"/"+JointCallingSparkConsts.GENOTYPE_GVCFS_PREFIX+this.idx);
        this.stateFile = this.jcOptions.output + "/" + JointCallingSparkConsts.JOB_STATE + "/" + JointCallingSparkConsts.GENOTYPE_GVCFS_PREFIX + this.idx + ".json";
    }

    public void submit() {
        if (isSuccess()) {
            return;
        }
        final String genotypePrefix = genotypeDir.getAbsolutePath()+"/"+JointCallingSparkConsts.GENOTYPE_GVCFS_PREFIX;
        // genotype partitions for each sub-interval, total partitions is about jcOptions.jobs*PARTITION_COEFFICIENT
        final int genotypePartitions = (int)Math.ceil(1.0*jcOptions.jobs*partitionCoeff/subIntervals.size());
        this.futures = new ArrayList<>();
        for (int j = 0; j < subIntervals.size(); ++j) {
            ArrayList<SimpleInterval> windows = SimpleIntervalUtils.splitIntervalByPartitionsAndBitSet(
                subIntervals.get(j), genotypePartitions, jcOptions.minVariantSites, subVariantBitSets.get(j));
            JavaRDD<SimpleInterval> windowsRDD = sc.parallelize(windows, windows.size());
            JavaFutureAction<List<String>> genotypeGVCFsFuture = windowsRDD
                .mapPartitionsWithIndex(new GVCFsSyncGenotyperSparkFunc(jcOptions.reference, combinedGVCFsList.get(j),
                    genotypeHeader, genotypePrefix+j+"-", jcOptions.dbsnp, jcOptions.genotypeArguments), false)
                .filter(x -> {return !x.equals("null");})
                .collectAsync();
            this.futures.add(genotypeGVCFsFuture);
        }
    }

    public List<String> get() {
        if (isSuccess()) {
            return load();
        }

        if (this.futures == null) {
            logger.error("Should run submit before get!");
            System.exit(1);
        }
        
        ArrayList<String> genotypeGVCFsList = new ArrayList<>();
        for (int j = 0; j < this.futures.size(); ++j) {
            try {
                genotypeGVCFsList.addAll(this.futures.get(j).get());
            } catch (Exception e) {
                logger.error("Failed to get genotyped gvcfs for interval: {}, {}", subIntervals.get(j), e.getMessage());
                System.exit(1);
            }
        }

        this.jobState.outPutFiles.put(0, genotypeGVCFsList);
        this.jobState.jobState = DPGTJobState.State.SUCCESS;

        if (jcOptions.deleteIntermediateResults) deleteFilesInRemoveList();
        writeStateFile();
        
        return genotypeGVCFsList;
    }

    public List<String> get(long timeout, TimeUnit unit) {
        if (isSuccess()) {
            return load();
        }

        if (this.futures == null) {
            logger.error("Should run submit before get!");
            System.exit(1);
        }
        
        ArrayList<String> genotypeGVCFsList = new ArrayList<>();
        for (int j = 0; j < this.futures.size(); ++j) {
            List<String> genotypeGVCF = null;
            try {
                genotypeGVCF = this.futures.get(j).get(timeout, unit);
            } catch (InterruptedException | ExecutionException | TimeoutException e) {
                if (TimeoutException.class.isInstance(e)) {
                    // time out
                    return null;
                } else {
                    logger.error("Failed to get genotyped gvcfs for interval: {}, {}", subIntervals.get(j), e.getMessage());
                    System.exit(1);
                }
            }
            if (genotypeGVCF != null) {
                genotypeGVCFsList.addAll(genotypeGVCF);
            }
        }

        this.jobState.outPutFiles.put(0, genotypeGVCFsList);
        this.jobState.jobState = DPGTJobState.State.SUCCESS;

        if (jcOptions.deleteIntermediateResults) deleteFilesInRemoveList();
        writeStateFile();
        
        return genotypeGVCFsList;
    }

    public List<String> load() {
        return this.jobState.outPutFiles.get(0);
    }

    public int getIndex() {
        return this.idx;
    }
}
