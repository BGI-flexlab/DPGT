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
#ifndef DPGT_GATK_VCF_CONSTANTS_HPP
#define DPGT_GATK_VCF_CONSTANTS_HPP

#include <string>
#include <vector>
#include "vcf_constants.hpp"

class GATKVCFConstants {
public:
    static const std::string CONTIG_ID_KEY;
    static const std::string CONTIG_LENGTH_KEY;
    static const std::string ASSEMBLY_NAME_KEY;

    //INFO keys
    static const std::string ALLELE_SPECIFIC_PREFIX;
    static const std::string AS_FILTER_STATUS_KEY;
    static const std::string RAW_RMS_MAPPING_QUALITY_DEPRECATED;  //NOTE: this is deprecated in favor of the new RAW_MQandDP below
    static const std::string MAPPING_QUALITY_DEPTH_DEPRECATED;  //NOTE: this is deprecated in favor of the new RAW_MQandDP below
    static const std::string RAW_MAPPING_QUALITY_WITH_DEPTH_KEY;
    static const std::string AS_RMS_MAPPING_QUALITY_KEY;
    static const std::string AS_RAW_RMS_MAPPING_QUALITY_KEY;
    static const std::string AS_CULPRIT_KEY;
    static const std::string AS_VQS_LOD_KEY;
    static const std::string ORIGINAL_AC_KEY; //SelectVariants
    static const std::string ORIGINAL_AF_KEY; //SelectVariants
    static const std::string ORIGINAL_AN_KEY; //SelectVariants
    static const std::string AC_ADJUSTED_KEY; //GnarlyGenotyper
    static const std::string BASE_QUAL_RANK_SUM_KEY;
    static const std::string BASE_QUAL_HISTOGRAM_KEY;
    static const std::string AS_BASE_QUAL_RANK_SUM_KEY;
    static const std::string AS_RAW_BASE_QUAL_RANK_SUM_KEY;
    static const std::string GENOTYPE_AND_VALIDATE_STATUS_KEY;
    static const std::string CLIPPING_RANK_SUM_KEY;
    static const std::string CULPRIT_KEY;
    static const std::string ORIGINAL_DP_KEY; //SelectVariants
    static const std::string DOWNSAMPLED_KEY;
    static const std::string EVENT_COUNT_IN_HAPLOTYPE_KEY; //M2
    static const std::string FISHER_STRAND_KEY;
    static const std::string AS_FISHER_STRAND_KEY;
    static const std::string AS_SB_TABLE_KEY;
    static const std::string SB_TABLE_KEY;
    static const std::string GQ_MEAN_KEY;
    static const std::string GQ_STDEV_KEY;
    static const std::string HAPLOTYPE_SCORE_KEY;
    static const std::string HI_CONF_DENOVO_KEY;
    static const std::string INTERVAL_GC_CONTENT_KEY;
    static const std::string INBREEDING_COEFFICIENT_KEY;
    static const std::string AS_INBREEDING_COEFFICIENT_KEY;
    static const std::string EXCESS_HET_KEY;
    static const std::string RAW_GENOTYPE_COUNT_KEY;
    static const std::string LIKELIHOOD_RANK_SUM_KEY;
    static const std::string LO_CONF_DENOVO_KEY;
    static const std::string MLE_ALLELE_COUNT_KEY;
    static const std::string MLE_ALLELE_FREQUENCY_KEY;
    static const std::string MAP_QUAL_RANK_SUM_KEY;
    static const std::string RAW_MAP_QUAL_RANK_SUM_KEY;
    static const std::string AS_MAP_QUAL_RANK_SUM_KEY;
    static const std::string AS_RAW_MAP_QUAL_RANK_SUM_KEY;
    static const std::string NOCALL_CHROM_KEY;
    static const std::string NUMBER_OF_DISCOVERED_ALLELES_KEY;
    static const std::string NEGATIVE_LABEL_KEY;
    static const std::string GENOTYPE_PRIOR_KEY;
    static const std::string POSITIVE_LABEL_KEY;
    static const std::string QUAL_BY_DEPTH_KEY;
    static const std::string AS_QUAL_BY_DEPTH_KEY;
    static const std::string AS_QUAL_KEY;
    static const std::string RAW_QUAL_APPROX_KEY;
    static const std::string AS_RAW_QUAL_APPROX_KEY;
    static const std::string VARIANT_DEPTH_KEY;
    static const std::string AS_VARIANT_DEPTH_KEY;
    static const std::string AS_ALT_ALLELE_DEPTH_KEY;
    static const std::string READ_POS_RANK_SUM_KEY;
    static const std::string AS_READ_POS_RANK_SUM_KEY;
    static const std::string AS_RAW_READ_POS_RANK_SUM_KEY;
    static const std::string REPEATS_PER_ALLELE_KEY;
    static const std::string REPEAT_UNIT_KEY;
    static const std::string SAMPLE_LIST_KEY;
    static const std::string STRAND_ODDS_RATIO_KEY;
    static const std::string AS_STRAND_ODDS_RATIO_KEY;
    static const std::string STR_PRESENT_KEY;
    static const std::string VQS_LOD_KEY;
    static const std::string CNN_1D_KEY;
    static const std::string CNN_2D_KEY;
    static const std::string F1R2_KEY;
    static const std::string F2R1_KEY;

    // Mutect2-specific INFO keys
    static const std::string TUMOR_LOG_10_ODDS_KEY;
    static const std::string NORMAL_LOG_10_ODDS_KEY;
    static const std::string IN_PON_KEY;
    static const std::string NORMAL_ARTIFACT_LOG_10_ODDS_KEY;
    static const std::string POPULATION_AF_KEY;
    static const std::string GERMLINE_QUAL_KEY;
    static const std::string SEQUENCING_QUAL_KEY;
    static const std::string POLYMERASE_SLIPPAGE_QUAL_KEY;
    static const std::string STRAND_QUAL_KEY;
    static const std::string CONTAMINATION_QUAL_KEY;
    static const std::string READ_ORIENTATION_QUAL_KEY;
    static const std::string ORIGINAL_CONTIG_MISMATCH_KEY;
    static const std::string N_COUNT_KEY;
    static const std::string AS_UNIQUE_ALT_READ_SET_COUNT_KEY;
    static const std::string MEDIAN_BASE_QUALITY_KEY;
    static const std::string MEDIAN_MAPPING_QUALITY_KEY;
    static const std::string MEDIAN_FRAGMENT_LENGTH_KEY;
    static const std::string MEDIAN_READ_POSITON_KEY;
    static const std::string UNITIG_SIZES_KEY;
    static const std::string ALIGNMENT_SCORE_DIFFERENCE_KEY;
    static const std::string JOINT_ALIGNMENT_COUNT_KEY;
    static const std::string REFERENCE_BASES_KEY;

    // Methylation-specific INFO Keys
    static const std::string UNCONVERTED_BASE_COVERAGE_KEY;
    static const std::string CONVERTED_BASE_COVERAGE_KEY;
    static const std::string METHYLATION_REFERENCE_CONTEXT_KEY;


    // FORMAT keys
    static const std::string ALLELE_BALANCE_KEY;
    static const std::string JOINT_LIKELIHOOD_TAG_NAME; //FamilyLikelihoodsUtils
    static const std::string JOINT_POSTERIOR_TAG_NAME; //FamilyLikelihoodsUtils
    const static std::string MIN_DP_FORMAT_KEY;
    static const std::string MAPPING_QUALITY_ZERO_BY_SAMPLE_KEY;
    static const std::string HAPLOTYPE_CALLER_PHASING_GT_KEY;
    static const std::string HAPLOTYPE_CALLER_PHASING_ID_KEY;
    static const std::string PHRED_SCALED_POSTERIORS_KEY; //FamilyLikelihoodsUtils / PosteriorLikelihoodsUtils
    static const std::string REFERENCE_GENOTYPE_QUALITY;
    static const std::string GENOTYPE_QUALITY_BY_ALLELE_BALANCE; //GnarlyGenotyper
    static const std::string GENOTYPE_QUALITY_BY_ALT_CONFIDENCE; //GnarlyGenotyper
    static const std::string STRAND_COUNT_BY_SAMPLE_KEY;
    static const std::string STRAND_BIAS_BY_SAMPLE_KEY;
    static const std::string FEATURIZED_READ_SETS_KEY;
    static const std::string HAPLOTYPE_EQUIVALENCE_COUNTS_KEY;
    static const std::string HAPLOTYPE_COMPLEXITY_KEY;
    static const std::string HAPLOTYPE_DOMINANCE_KEY;
    const static std::string TRANSMISSION_PROBABILITY_KEY; //PhaseByTransmission
    static const std::string FRAGMENT_ALLELE_DEPTHS;

    // M2-specific FORMAT keys
    static const std::string ALLELE_FRACTION_KEY;

    //FILTERS
    /* Note that many filters used throughout GATK (most notably in VariantRecalibration) are dynamic,
       their names (or descriptions) depend on some threshold.  Those filters are not included here
     */
    static const std::string CLUSTERED_EVENTS_FILTER_NAME; //M2
    static const std::string GERMLINE_RISK_FILTER_NAME; //M2
    static const std::string LOW_QUAL_FILTER_NAME;
    static const std::string ALIGNMENT_ARTIFACT_FILTER_NAME;
    static const std::string PON_FILTER_NAME; //M2
    static const std::string POLYMERASE_SLIPPAGE; //M2
    static const std::string TUMOR_EVIDENCE_FILTER_NAME; //M2
    static const std::string MULTIALLELIC_FILTER_NAME; //M2
    static const std::string STRAND_ARTIFACT_FILTER_NAME; // M2
    static const std::string DUPLICATED_EVIDENCE_FILTER_NAME;
    const static std::string ARTIFACT_IN_NORMAL_FILTER_NAME;
    const static std::string MEDIAN_BASE_QUALITY_FILTER_NAME;
    const static std::string MEDIAN_MAPPING_QUALITY_FILTER_NAME;
    const static std::string MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME;
    const static std::string READ_POSITION_FILTER_NAME;
    const static std::string CONTAMINATION_FILTER_NAME;
    const static std::string READ_ORIENTATION_ARTIFACT_FILTER_NAME;
    const static std::string BAD_HAPLOTYPE_FILTER_NAME;
    const static std::string STRICT_STRAND_BIAS_FILTER_NAME;
    const static std::string N_RATIO_FILTER_NAME;
    const static std::string ALLELE_FRACTION_FILTER_NAME;
    static const std::string POSSIBLE_NUMT_FILTER_NAME;
    static const std::string LOW_HET_FILTER_NAME;
    static const std::string FAIL;
    static const std::string SITE_LEVEL_FILTERS;


    static const std::vector<std::string> MUTECT_FILTER_NAMES;

    static const std::vector<std::string> MUTECT_AS_FILTER_NAMES;

    // Symbolic alleles
    const static std::string SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG;
    const static std::string NON_REF_SYMBOLIC_ALLELE_NAME;
    static const std::string ALLELE_SPECIFIC_ANNOTATION_PREFIX;
    static const std::string PSEUDO_DEPTH_KEY;
    static const std::string PSEUDO_FRACTION_KEY;

    // INFO lines
    static const std::string RAW_MAPPING_QUALITY_WITH_DEPTH_KEY_LINE;
};



#endif  // DPGT_GATK_VCF_CONSTANTS_HPP
