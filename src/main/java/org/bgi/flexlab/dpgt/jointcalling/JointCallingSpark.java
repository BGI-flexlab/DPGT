package org.bgi.flexlab.dpgt.jointcalling;

import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.bgi.flexlab.dpgt.utils.SimpleIntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.SparkConf;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.io.InputStream;

import org.apache.log4j.PropertyConfigurator;
import org.apache.commons.io.FileUtils;

public class JointCallingSpark {
    private static final Logger logger = LoggerFactory.getLogger(JointCallingSpark.class);

    private static final String HEADER_DIR = "header";  // outdir/header
    private static final String GENOTYPE_HEADER = "genotype_header.vcf.gz";  // outdir/header/genotype_header.vcf.gz
    private static final String COMBINE_GVCFS_PREFIX = "combine."; // outdir/combine.idx
    private static final String GENOTYPE_GVCFS_PREFIX = "genotype."; // outdir/genotype.idx
    private static final String OUTPUT_NAME = "result.vcf.gz"; // outdir/result.vcf.gz

    private static final int GENOTYPE_PARTITION_COEFFICIENT = 4;

    static{
        try{
            System.loadLibrary("cdpgt");
        }catch(UnsatisfiedLinkError e){
            logger.error("Failed to load Native code library cdpgt.\n" + e);
            System.exit(1);
        }
    }

    public static void main(String[] args) throws Exception {
        ClassLoader classloader = Thread.currentThread().getContextClassLoader();
        InputStream log4jPropertiestIs = classloader.getResourceAsStream("log4j.properties");
        PropertyConfigurator.configure(log4jPropertiestIs);
        JointCallingSparkOptions jcOptions = new JointCallingSparkOptions();
        jcOptions.parse(args);
        List<SimpleInterval> intervalsToTravers = jcOptions.getIntervalsToTravers();

        SparkConf conf = new SparkConf().setAppName("VariantSiteFinder");
        JavaSparkContext sc = new JavaSparkContext(conf);

        JavaRDD<String> vcfpathsRDDPartitionByCombineParts = sc.textFile(jcOptions.input, jcOptions.numCombinePartitions);

        File headerDir = new File(jcOptions.output+"/"+HEADER_DIR);
        if (!headerDir.exists()) {
            headerDir.mkdirs();
        }
        combineVCFHeaders(vcfpathsRDDPartitionByCombineParts, headerDir.getAbsolutePath(), jcOptions.genotypeArguments);
        final String genotypeHeader = headerDir.getAbsolutePath() + "/" + GENOTYPE_HEADER;

        JavaRDD<String> vcfpathsRDDPartitionByJobs = sc.textFile(jcOptions.input, jcOptions.jobs);

        final String outputPath = jcOptions.output + "/" + OUTPUT_NAME;

        for (int i = 0; i < intervalsToTravers.size(); ++i) {
            SimpleInterval interval = intervalsToTravers.get(i);
            logger.info("Cycle {}/{}", i+1, intervalsToTravers.size());
            logger.info("Processing interval: {}", interval.toString());
            logger.info("Finding variant sites in {}", interval.toString());
            BitSet variantSiteSetData = vcfpathsRDDPartitionByJobs.
                mapPartitionsWithIndex(new VariantSiteFinderSparkFunc(interval.getContig(), interval.getStart()-1, interval.getEnd()-1), false).
                reduce((x, y) -> {x.or(y); return x;});
            if (variantSiteSetData.isEmpty()) {
                logger.info("skip interval: {}, because there is no variant site in it.", interval.toString());
                continue;
            }

            Broadcast<byte[]> variantSiteSetBytesBc = sc.broadcast(variantSiteSetData.toByteArray());

            logger.info("Combining gvcfs in {}", interval.toString());
            File combineDir = new File(jcOptions.output+"/"+COMBINE_GVCFS_PREFIX+i);
            if (!combineDir.exists()) {
                combineDir.mkdirs();
            }
            List<String> combinedGVCFs = vcfpathsRDDPartitionByCombineParts.mapPartitionsWithIndex(new CombineGVCFsOnSitesSparkFunc(
                jcOptions.reference, combineDir.getAbsolutePath()+"/"+COMBINE_GVCFS_PREFIX, interval, variantSiteSetBytesBc), false)
                .filter(x -> {return !x.equals("null");})
                .collect();
            
            variantSiteSetBytesBc.unpersist();  // Asynchronously delete cached copies of this broadcast on the executors.

            // List<String> combinedGVCFs = Arrays.asList(combineDir.list((x, y) ->{return y.endsWith("vcf.gz");}));
            // for (int j = 0; j < combinedGVCFs.size(); ++j) {
            //     combinedGVCFs.set(j, combineDir.getAbsolutePath() + "/" + combinedGVCFs.get(j));
            // }

            logger.info("Genotyping gvcfs in {}", interval.toString());
            File genotypeDir = new File(jcOptions.output+"/"+GENOTYPE_GVCFS_PREFIX+i);
            if (!genotypeDir.exists()) {
                genotypeDir.mkdirs();
            }
            final String genotypePrefix = genotypeDir.getAbsolutePath()+"/"+GENOTYPE_GVCFS_PREFIX;
            ArrayList<SimpleInterval> windows = SimpleIntervalUtils.splitIntervalByPartitions(interval, jcOptions.jobs * GENOTYPE_PARTITION_COEFFICIENT);
            JavaRDD<SimpleInterval> windowsRDD = sc.parallelize(windows, windows.size());
            List<String> genotypeGVCFs = windowsRDD.
                mapPartitionsWithIndex(new GVCFsSyncGenotyperSparkFunc(jcOptions.reference, combinedGVCFs,
                    genotypeHeader, genotypePrefix, jcOptions.dbsnp, jcOptions.genotypeArguments), false)
                .filter(x -> {return !x.equals("null");})
                .collect();
            
            logger.info("Concating genotyped vcfs in {} to result", interval.toString());
            JavaRDD<String> concateGVCFsRDD = sc.parallelize(genotypeGVCFs, 1);
            concateGVCFsRDD.mapPartitionsWithIndex(new ConcatGenotypeGVCFsSparkFunc((i==0 ? genotypeHeader : null), outputPath), false).collect();

            if (jcOptions.deleteIntermediateResults) {
                FileUtils.deleteDirectory(combineDir);
                FileUtils.deleteDirectory(genotypeDir);
            }
        }

        if (jcOptions.deleteIntermediateResults) {
            FileUtils.deleteDirectory(headerDir);
        }

        sc.close();
    }

    private static VCFHeader combineVCFHeaders(JavaRDD<String> vcfpathsRDD, final String outdir, final GenotypeCalculationArgumentCollection genotypeArgs) {
        List<String> combinedHeaders = vcfpathsRDD.mapPartitionsWithIndex(new CombineVCFHeadersSparkFunc(outdir), false)
            .filter(x -> {return !x.equals("null");})
            .collect();

        // combine headers of each partition to generate combined header for all input vcfs
        VCFHeaderCombiner combiner = new VCFHeaderCombiner();
        String combinedHeader = outdir+"/header.vcf.gz";
        String[] combinedHeadersArray = new String[combinedHeaders.size()];
        combinedHeaders.toArray(combinedHeadersArray);
        combiner.Combine(combinedHeadersArray, combinedHeader);

        // make combined header for genotyping by add genotyping specific headers
        GVCFsSyncGenotyper genotyper = new GVCFsSyncGenotyper();
        VCFHeader combinedHeaderForGenotyping = genotyper.makeCombinedHeaderForGenotyping(combinedHeader, genotypeArgs);
        String combinedHeaderForGenotypingPath = outdir+"/"+GENOTYPE_HEADER;
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        builder.setOutputFile(combinedHeaderForGenotypingPath);
        builder.setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF);
        VariantContextWriter writer = builder.build();
        writer.writeHeader(combinedHeaderForGenotyping);
        writer.close();
        
        return combinedHeaderForGenotyping;
    }
}
