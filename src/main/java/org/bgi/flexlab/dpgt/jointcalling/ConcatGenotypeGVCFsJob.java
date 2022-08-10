package org.bgi.flexlab.dpgt.jointcalling;

import java.util.List;
import java.io.File;
import java.util.Arrays;
import org.bgi.flexlab.dpgt.utils.DPGTJob;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;


public class ConcatGenotypeGVCFsJob extends DPGTJob<Integer> {
    private final JointCallingSparkOptions jcOptions;
    private JavaSparkContext sc;
    private List<String> genotypeGVCFsList;
    private String genotypeHeader;
    private String outputPath;
    public ConcatGenotypeGVCFsJob(final JointCallingSparkOptions jcOptions, final JavaSparkContext sc, final List<String> genotypeGVCFsList, final String genotypeHeader, final String outputPath) {
        this.jcOptions = jcOptions;
        this.sc = sc;
        this.genotypeGVCFsList = genotypeGVCFsList;
        this.genotypeHeader = genotypeHeader;
        File outputFile = new File(outputPath);
        this.outputPath = outputFile.getAbsolutePath();
        this.stateFile = this.jcOptions.output + "/" + JointCallingSparkConsts.JOB_STATE + "/" + JointCallingSparkConsts.CONCAT_GENOTYPE_GVCFS;
    }

    public ConcatGenotypeGVCFsJob(final JointCallingSparkOptions jcOptions) {
        this.jcOptions = jcOptions;
        this.stateFile = this.jcOptions.output + "/" + JointCallingSparkConsts.JOB_STATE + "/" + JointCallingSparkConsts.CONCAT_GENOTYPE_GVCFS;
    }

    public Integer work() {
        JavaRDD<String> concateGVCFsRDD = sc.parallelize(genotypeGVCFsList, 1);
        concateGVCFsRDD.mapPartitionsWithIndex(new ConcatGenotypeGVCFsSparkFunc(genotypeHeader, outputPath), false).collect();
        this.jobState.outPutFiles.put(0, Arrays.asList(new String[]{outputPath}));
        return 0;
    }

    public Integer load() {
        return 0;
    }
}
