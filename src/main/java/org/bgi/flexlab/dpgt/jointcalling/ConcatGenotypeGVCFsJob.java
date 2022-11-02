package org.bgi.flexlab.dpgt.jointcalling;

import java.util.ArrayList;
import java.util.List;
import java.io.File;
import java.io.FileOutputStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.bgi.flexlab.dpgt.utils.DPGTJobAsync;
import org.bgi.flexlab.dpgt.utils.DPGTJobState;
import org.apache.spark.api.java.JavaFutureAction;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;


public class ConcatGenotypeGVCFsJob extends DPGTJobAsync<List<String>, Integer> {
    private static final Logger logger = LoggerFactory.getLogger(ConcatGenotypeGVCFsJob.class);
    private final JointCallingSparkOptions jcOptions;
    private JavaSparkContext sc;
    private List<String> genotypeGVCFsList;
    private String genotypeHeader;
    private String outputPath;
    private int idx;
    public ConcatGenotypeGVCFsJob(final JointCallingSparkOptions jcOptions, final JavaSparkContext sc,
        final List<String> genotypeGVCFsList, final String genotypeHeader, final String outputPath, int idx)
    {
        this.jcOptions = jcOptions;
        this.sc = sc;
        this.genotypeGVCFsList = genotypeGVCFsList;
        this.genotypeHeader = genotypeHeader;
        File outputFile = new File(outputPath);
        this.outputPath = outputFile.getAbsolutePath();
        this.idx = idx;
        this.stateFile = String.format("%s/%s/%s-%d.json",
            this.jcOptions.output, JointCallingSparkConsts.JOB_STATE,
            JointCallingSparkConsts.CONCAT_GENOTYPE_STATE_FILE_PREFIX, this.idx);
    }

    public ConcatGenotypeGVCFsJob(final JointCallingSparkOptions jcOptions, int idx) {
        this.jcOptions = jcOptions;
        this.idx = idx;
        this.stateFile = String.format("%s/%s/%s-%d.json",
            this.jcOptions.output, JointCallingSparkConsts.JOB_STATE,
            JointCallingSparkConsts.CONCAT_GENOTYPE_STATE_FILE_PREFIX, this.idx);
    }

    public void submit() {
        if (isSuccess()) {
            return;
        }
        JavaRDD<String> concateGVCFsRDD = sc.parallelize(genotypeGVCFsList, 1);
        JavaFutureAction<List<String>> futureAction = concateGVCFsRDD.mapPartitionsWithIndex(
            new ConcatGenotypeGVCFsSparkFunc(genotypeHeader, outputPath), false).collectAsync();
        this.futures = new ArrayList<>();
        this.futures.add(futureAction);
    }

    public Integer get() {
        if (isSuccess()) {
            return load();
        }

        if (this.futures == null) {
            logger.error("Should run submit before get!");
            System.exit(1);
        }

        for (int i = 0; i < this.futures.size(); ++i) {
            try {
                this.futures.get(i).get();
            } catch (Exception e) {
                logger.error("Failed to get concatenate vcf, {}", e.getMessage());
                System.exit(1);
            }
        }

        // create an empty flag file
        File jobSuccessFlag = new File(String.format("%s.%s_%d", outputPath, JointCallingSparkConsts.JOB_SUCCESS_FLAG, this.idx));
        try {
            FileOutputStream jobSuccessFlagStream = new FileOutputStream(jobSuccessFlag);
            jobSuccessFlagStream.close();
        } catch (Exception e) {
            logger.error("Failed to get concatenate vcf, {}", e.getMessage());
            System.exit(1);
        }

        ArrayList<String> outPutFiles = new ArrayList<>();
        outPutFiles.add(outputPath);
        outPutFiles.add(jobSuccessFlag.getAbsolutePath());

        this.jobState.outPutFiles.put(0, outPutFiles);
        this.jobState.jobState = DPGTJobState.State.SUCCESS;

        return 0;
    }

    public Integer load() {
        return 0;
    }
}
