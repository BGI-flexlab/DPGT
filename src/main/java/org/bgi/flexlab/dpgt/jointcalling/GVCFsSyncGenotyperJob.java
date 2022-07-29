package org.bgi.flexlab.dpgt.jointcalling;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.bgi.flexlab.dpgt.utils.DPGTJob;
import org.bgi.flexlab.dpgt.utils.SimpleIntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.JavaFutureAction;

public class GVCFsSyncGenotyperJob extends DPGTJob<List<String>> {
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
        this.stateFile = this.jcOptions.output + "/" + JointCallingSparkConsts.JOB_STATE + "/" + JointCallingSparkConsts.GENOTYPE_GVCFS_PREFIX + this.idx + ".json";
    }

    public List<String> work() {
        final String genotypePrefix = genotypeDir.getAbsolutePath()+"/"+JointCallingSparkConsts.GENOTYPE_GVCFS_PREFIX;
        // genotype partitions for each sub-interval, total partitions is about jcOptions.jobs*PARTITION_COEFFICIENT
        final int genotypePartitions = (int)Math.ceil(1.0*jcOptions.jobs*partitionCoeff/subIntervals.size());
        ArrayList<JavaFutureAction<List<String>>> genotypeGVCFsFutures = new ArrayList<>();
        for (int j = 0; j < subIntervals.size(); ++j) {
            ArrayList<SimpleInterval> windows = SimpleIntervalUtils.splitIntervalByPartitionsAndBitSet(
                subIntervals.get(j), genotypePartitions, subVariantBitSets.get(j));
            JavaRDD<SimpleInterval> windowsRDD = sc.parallelize(windows, windows.size());
            JavaFutureAction<List<String>> genotypeGVCFsFuture = windowsRDD
                .mapPartitionsWithIndex(new GVCFsSyncGenotyperSparkFunc(jcOptions.reference, combinedGVCFsList.get(j),
                    genotypeHeader, genotypePrefix+j+"-", jcOptions.dbsnp, jcOptions.genotypeArguments), false)
                .filter(x -> {return !x.equals("null");})
                .collectAsync();
            genotypeGVCFsFutures.add(genotypeGVCFsFuture);
        }
        
        ArrayList<String> genotypeGVCFsList = new ArrayList<>();
        for (int j = 0; j < genotypeGVCFsFutures.size(); ++j) {
            try {
                genotypeGVCFsList.addAll(genotypeGVCFsFutures.get(j).get());
            } catch (Exception e) {
                logger.error("Failed to get genotyped gvcfs for interval: {}, {}", subIntervals.get(j), e.getMessage());
                System.exit(1);
            }
        }

        this.jobState.outPutFiles.put(0, genotypeGVCFsList);

        return genotypeGVCFsList;
    }

    public List<String> load() {
        return this.jobState.outPutFiles.get(0);
    }
}
