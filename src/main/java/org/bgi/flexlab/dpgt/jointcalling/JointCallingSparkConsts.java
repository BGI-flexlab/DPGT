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
    public static final String JOB_STATE = "job_state";                     // outdir/job_state
    public static final Integer GET_FUTURE_TIMEOUT = 100;                   // max time out value in milliseconds
    public static final String COMBINE_HEADER_STATE_FILE = "combine_vcf_header.json";      // outdir/job_state/combine_vcf_header.json
    public static final String CONCAT_GENOTYPE_STATE_FILE_PREFIX = "concat_genotype_gvcfs";       // outdir/job_state/concat_genotype_gvcfs.json
    public static final String JOB_SUCCESS_FLAG = "SUCCESS";
    public static final String RESULT_VCF_FILE_SIZE_KEY = "result_vcf_file_size";   // vcf file size after one concat vcf job done
}
