package org.bgi.flexlab.dpgt.jointcalling;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.io.File;
import java.io.FileOutputStream;
import java.nio.file.Paths;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.bgi.flexlab.dpgt.utils.DPGTJobAsync;
import org.bgi.flexlab.dpgt.utils.DPGTJobState;
import org.apache.spark.api.java.JavaFutureAction;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;


public class ConcatGenotypeGVCFsJob extends DPGTJobAsync<List<String>, Integer> {
    private static final Logger logger = LoggerFactory.getLogger(ConcatGenotypeGVCFsJob.class);
    private boolean isNull = false;
    private JointCallingSparkOptions jcOptions;
    private JavaSparkContext sc;
    private List<String> genotypeGVCFsList;
    private String genotypeHeader;
    private String outputPath;
    private int idx;
    public ConcatGenotypeGVCFsJob(final JointCallingSparkOptions jcOptions, final JavaSparkContext sc,
        final List<String> genotypeGVCFsList, final String genotypeHeader, int idx)
    {
        this.jcOptions = jcOptions;
        this.sc = sc;
        this.genotypeGVCFsList = genotypeGVCFsList;
        this.genotypeHeader = genotypeHeader;
        File outputFile = new File(this.jcOptions.outputPrefix + "." + Integer.toString(idx) + "." + JointCallingSparkConsts.OUTPUT_SUFFIX);
        this.outputPath = outputFile.getAbsolutePath();
        this.idx = idx;
        this.stateFile = getStateFilePath();
    }

    public ConcatGenotypeGVCFsJob(final JointCallingSparkOptions jcOptions, int idx) {
        this.jcOptions = jcOptions;
        this.idx = idx;
        this.stateFile = getStateFilePath();
    }

    public ConcatGenotypeGVCFsJob() {
        this.isNull = true;
    }

    boolean isNull() {
        return this.isNull;
    }
    
    void setToNull() {
        this.isNull = true;
        this.jcOptions = null;
        this.sc = null;
        this.genotypeGVCFsList = null;
        this.genotypeHeader = null;
        this.outputPath = null;
        this.idx = 0;
    }

    public void submit() {
        if (isSuccess()) {
            return;
        }
        JavaRDD<String> concateGVCFsRDD = sc.parallelize(genotypeGVCFsList, 1);
        JavaFutureAction<List<String>> futureAction = concateGVCFsRDD.mapPartitionsWithIndex(
            new ConcatGenotypeGVCFsSparkFunc(genotypeHeader, outputPath, getOutputFileSizeBefore()), false).collectAsync();
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

        File outputFile = new File(outputPath);
        this.jobState.metaData.put(JointCallingSparkConsts.RESULT_VCF_FILE_SIZE_KEY, String.valueOf(outputFile.length()));

        this.jobState.jobState = DPGTJobState.State.SUCCESS;

        if (jcOptions.deleteIntermediateResults) deleteFilesInRemoveList();
        writeStateFile();

        return 0;
    }

    public Integer get(long timeout, TimeUnit unit) {
        if (isSuccess()) {
            return load();
        }

        if (this.futures == null) {
            logger.error("Should run submit before get!");
            System.exit(1);
        }

        for (int i = 0; i < this.futures.size(); ++i) {
            try {
                this.futures.get(i).get(timeout, unit);
            } catch (InterruptedException | ExecutionException | TimeoutException e) {
                if (TimeoutException.class.isInstance(e)) {
                    return null;
                } else {
                    logger.error("Failed to get concatenate vcf, {}", e.getMessage());
                    System.exit(1);
                }
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

        File outputFile = new File(outputPath);
        this.jobState.metaData.put(JointCallingSparkConsts.RESULT_VCF_FILE_SIZE_KEY, String.valueOf(outputFile.length()));

        this.jobState.jobState = DPGTJobState.State.SUCCESS;

        if (jcOptions.deleteIntermediateResults) deleteFilesInRemoveList();
        writeStateFile();

        return 0;
    }

    public Integer load() {
        return 0;
    }

    /**
     * get output result vcf file size before this job
     * @return
     */
    private long getOutputFileSizeBefore() {
        if (this.idx == 0) {
            return 0;
        } else {
            String preStateFile = String.format("%s/%s/%s.%d.json",
                this.jcOptions.output, JointCallingSparkConsts.JOB_STATE,
                JointCallingSparkConsts.CONCAT_GENOTYPE_STATE_FILE_PREFIX, this.idx-1);
            DPGTJobState preJobState = readStateFile(preStateFile);
            String fileSizeStr = preJobState.metaData.get(JointCallingSparkConsts.RESULT_VCF_FILE_SIZE_KEY);
            if (fileSizeStr == null) {
                logger.error("key {} not in state file {} metadata", JointCallingSparkConsts.RESULT_VCF_FILE_SIZE_KEY, preStateFile);
                System.exit(1);
            }
            long fileSize = 0;
            try {
                fileSize = Long.parseLong(fileSizeStr);
                return fileSize;
            } catch (NumberFormatException e) {
                logger.error("{} state file metadata key {} not have a string value that can be convert to long",
                    JointCallingSparkConsts.RESULT_VCF_FILE_SIZE_KEY, preStateFile);
                System.exit(1);
            }
            return fileSize;
        }
    }

    private String getStateFilePath() {
        return String.format("%s/%s/%s.%d.json",
            this.jcOptions.output, JointCallingSparkConsts.JOB_STATE,
            JointCallingSparkConsts.CONCAT_GENOTYPE_STATE_FILE_PREFIX, this.idx);
    }
}
