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

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.io.File;
import java.io.FileOutputStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.bgi.flexlab.dpgt.utils.DPGTJobAsync;
import org.bgi.flexlab.dpgt.utils.DPGTJobState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
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
    private int idx;
    private SimpleInterval interval;
    private String outputPath;
    public ConcatGenotypeGVCFsJob(final JointCallingSparkOptions jcOptions, final JavaSparkContext sc,
        final List<String> genotypeGVCFsList, final String genotypeHeader, int idx, SimpleInterval interval)
    {
        this.jcOptions = jcOptions;
        this.sc = sc;
        this.genotypeGVCFsList = genotypeGVCFsList;
        this.genotypeHeader = genotypeHeader;
        this.idx = idx;
        this.interval = interval;
        String outputFileName = String.format("%s.%s.%d.%s", this.jcOptions.outputPrefix,
            getSimpleIntervalString(), this.idx, JointCallingSparkConsts.OUTPUT_SUFFIX);
        File outputFile = new File(outputFileName);
        this.outputPath = outputFile.getAbsolutePath();
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
        return 0;
    }

    private String getStateFilePath() {
        return String.format("%s/%s/%s.%d.json",
            this.jcOptions.output, JointCallingSparkConsts.JOB_STATE,
            JointCallingSparkConsts.CONCAT_GENOTYPE_STATE_FILE_PREFIX, this.idx);
    }

    private String getSimpleIntervalString() {
        int seqLen = jcOptions.getSequenceDict().getSequence(this.interval.getContig()).getSequenceLength();
        // if interval covers the whole chromosome, return the contig name
        if (this.interval.getStart() == 1 && this.interval.getEnd() == seqLen) {
            return this.interval.getContig();
        } else {
            return String.format("%s_%d_%d", interval.getContig(), interval.getStart(), interval.getEnd());
        }
    }
}
