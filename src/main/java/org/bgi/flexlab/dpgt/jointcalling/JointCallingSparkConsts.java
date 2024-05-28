package org.bgi.flexlab.dpgt.jointcalling;


public class JointCallingSparkConsts {

    public static final String HEADER_DIR = "header";                       // outdir/header
    public static final String COMBINE_HEADER = "header.vcf.gz";            // outdir/header/header.vcf.gz
    public static final String GENOTYPE_HEADER = "genotype_header.vcf.gz";  // outdir/header/genotype_header.vcf.gz
    public static final String VARIANT_SITE_PREFIX = "variant.";            // outdir/variant.idx.site
    public static final String VARIANT_SITE_SUFFIX = ".site";               // outdir/variant.idx.site
    public static final String COMBINE_GVCFS_PREFIX = "combine.";           // outdir/combine.idx
    public static final String GENOTYPE_GVCFS_PREFIX = "genotype.";         // outdir/genotype.idx
    public static final String OUTPUT_NAME = "result.vcf.gz";               // outdir/result.vcf.gz
    public static final String OUTPUT_PREFIX = "result";                    // outdir/result.*.vcf.gz
    public static final String OUTPUT_SUFFIX = "vcf.gz";                    // outdir/result.*.vcf.gz
    public static final String VCF_PAIRS = "vcf_and_index_pairs.list";      // outdir/vcf_and_index_pairs.list
    public static final String JOB_STATE = "job_state";                     // outdir/job_state
    public static final Integer GET_FUTURE_TIMEOUT = 100;                   // max time out value in milliseconds
    public static final String COMBINE_HEADER_STATE_FILE = "combine_vcf_header.json";      // outdir/job_state/combine_vcf_header.json
    public static final String CONCAT_GENOTYPE_STATE_FILE_PREFIX = "concat_genotype_gvcfs";       // outdir/job_state/concat_genotype_gvcfs.json
    public static final String JOB_SUCCESS_FLAG = "SUCCESS";
    public static final String RESULT_VCF_FILE_SIZE_KEY = "result_vcf_file_size";   // vcf file size after one concat vcf job done
}
