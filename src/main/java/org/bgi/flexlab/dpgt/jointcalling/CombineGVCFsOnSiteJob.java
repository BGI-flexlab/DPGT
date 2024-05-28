package org.bgi.flexlab.dpgt.jointcalling;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.cloud.hadoop.gcsio.testing.InMemoryBucketEntry;

import org.bgi.flexlab.dpgt.utils.DPGTJob;
import org.bgi.flexlab.dpgt.utils.VariantSiteSetUtils;
import org.bgi.flexlab.dpgt.utils.SimpleIntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.api.java.JavaFutureAction;

public class CombineGVCFsOnSiteJob extends DPGTJob<Integer> {
    private static final Logger logger = LoggerFactory.getLogger(CombineGVCFsOnSiteJob.class);
    private JavaRDD<String> vcfPairsRDDPartitionByCombineParts;
    private JointCallingSparkOptions jcOptions;
    private JavaSparkContext sc;
    private int idx;
    private SimpleInterval interval;
    private int partitionCoeff;
    private BitSet variantSiteSetData;
    public File combineDir;
    // results
    public ArrayList<SimpleInterval> subIntervals = null;
    public ArrayList<BitSet> subVariantBitSets = null;
    public ArrayList<List<String>> combinedGVCFsList = null;

    public CombineGVCFsOnSiteJob(final JavaRDD<String> vcfPairsRDDPartitionByCombineParts, final JointCallingSparkOptions jcOptions,
        final JavaSparkContext sc, final int idx, final SimpleInterval interval, final int partitionCoeff, final BitSet variantSiteSetData)
    {
        this.vcfPairsRDDPartitionByCombineParts = vcfPairsRDDPartitionByCombineParts;
        this.jcOptions = jcOptions;
        this.sc = sc;
        this.idx = idx;
        this.interval = interval;
        this.combineDir = new File(jcOptions.output+"/"+JointCallingSparkConsts.COMBINE_GVCFS_PREFIX+this.idx);
        if (!this.combineDir.exists()) {
            this.combineDir.mkdirs();
        }
        this.partitionCoeff = this.jcOptions.useLix ? 1 : partitionCoeff;
        this.variantSiteSetData = variantSiteSetData;
        this.stateFile = this.jcOptions.output + "/" + JointCallingSparkConsts.JOB_STATE + "/" + JointCallingSparkConsts.COMBINE_GVCFS_PREFIX + this.idx + ".json";
    } 

    public Integer work() {
        if (this.jcOptions.useLix) {
            return workNoSubintervals();
        } else {
            return workWithSubIntervals();
        }
    }

    public Integer load() {
        if (this.jcOptions.useLix) {
            // not split interval into sub intervals
            this.subIntervals = new ArrayList<>();
            this.subIntervals.add(interval);
            this.subVariantBitSets = new ArrayList<>();
            this.subVariantBitSets.add(variantSiteSetData);
        } else {
            // split interval into sub intervals
            this.subIntervals = SimpleIntervalUtils.splitIntervalByPartitionsAndBitSet(interval,
                Math.max((int)Math.ceil(partitionCoeff*1.0*jcOptions.jobs/jcOptions.numCombinePartitions), 1),
                jcOptions.minVariantSites, variantSiteSetData);
            this.subVariantBitSets = VariantSiteSetUtils.getSubBitSets(variantSiteSetData, interval, subIntervals);
            this.combinedGVCFsList = new ArrayList<>();
        }

        for (final List<String> value: this.jobState.outPutFiles.values()) {
            this.combinedGVCFsList.add(value);
        }

        return 0;
    }

    private Integer workNoSubintervals() {
        // not split interval into sub intervals
        this.subIntervals = new ArrayList<>();
        this.subIntervals.add(interval);
        this.subVariantBitSets = new ArrayList<>();
        this.subVariantBitSets.add(variantSiteSetData);
        return workCommon();
    }

    private Integer workWithSubIntervals() {
        this.subIntervals = SimpleIntervalUtils.splitIntervalByPartitionsAndBitSet(interval,
            Math.max((int)Math.ceil(partitionCoeff*1.0*jcOptions.jobs/jcOptions.numCombinePartitions), 1),
            jcOptions.minVariantSites, variantSiteSetData);
        this.subVariantBitSets = VariantSiteSetUtils.getSubBitSets(variantSiteSetData, interval, subIntervals);
        return workCommon();
    }

    private Integer workCommon() {
        ArrayList<Broadcast<byte[]>> subVariantBitSetsBcs = new ArrayList<>();
        for (BitSet b: subVariantBitSets)
        {
            subVariantBitSetsBcs.add(sc.broadcast(b.toByteArray()));
        }

        final String combinePrefix = combineDir.getAbsolutePath()+"/"+JointCallingSparkConsts.COMBINE_GVCFS_PREFIX;
        ArrayList<JavaFutureAction<List<String>>> combinedGVCFsFutures = new ArrayList<>();
        for (int j = 0; j < subIntervals.size(); ++j) {
            JavaFutureAction<List<String>> combinedGVCFsFuture = vcfPairsRDDPartitionByCombineParts.mapPartitionsWithIndex(new CombineGVCFsOnSitesSparkFunc(
                jcOptions.reference, combinePrefix+j+"-", this.jcOptions.useLix ? 1 : 0, subIntervals.get(j), subVariantBitSetsBcs.get(j)), false)
                .filter(x -> {return !x.equals("null");})
                .collectAsync();  // collect async to run on multiple sub intervals at the same time
            combinedGVCFsFutures.add(combinedGVCFsFuture);
        }
        this.combinedGVCFsList = new ArrayList<>();
        for (int i = 0; i < combinedGVCFsFutures.size(); ++i) {
            try {
                combinedGVCFsList.add(combinedGVCFsFutures.get(i).get());
            } catch (Exception e) {
                logger.error("Failed to get combined gvcfs for interval: {}, {}", subIntervals.get(i), e.getMessage());
                System.exit(1);
            }
        }

        for (Broadcast<byte[]> bc: subVariantBitSetsBcs) {
            bc.unpersist();
        }

        for (int i = 0; i < combinedGVCFsList.size(); ++i) {
            this.jobState.outPutFiles.put(i, new ArrayList<>());
            this.jobState.outPutFiles.get(i).addAll(combinedGVCFsList.get(i));
        }

        return 0;
    }
}
