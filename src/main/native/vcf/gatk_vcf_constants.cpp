#include "gatk_vcf_constants.hpp"
#include <vector>


const std::string GATKVCFConstants::CONTIG_ID_KEY =                      "ID";
const std::string GATKVCFConstants::CONTIG_LENGTH_KEY =                  "length";
const std::string GATKVCFConstants::ASSEMBLY_NAME_KEY =                  "assembly";

//INFO keys
const std::string GATKVCFConstants::ALLELE_SPECIFIC_PREFIX =             "AS_";
const std::string GATKVCFConstants::AS_FILTER_STATUS_KEY =               "AS_FilterStatus";
const std::string GATKVCFConstants::RAW_RMS_MAPPING_QUALITY_DEPRECATED =        "RAW_MQ";  //NOTE: this is deprecated in favor of the new RAW_MQandDP below
const std::string GATKVCFConstants::MAPPING_QUALITY_DEPTH_DEPRECATED =   "MQ_DP";  //NOTE: this is deprecated in favor of the new RAW_MQandDP below
const std::string GATKVCFConstants::RAW_MAPPING_QUALITY_WITH_DEPTH_KEY = "RAW_MQandDP";
const std::string GATKVCFConstants::AS_RMS_MAPPING_QUALITY_KEY =         "AS_MQ";
const std::string GATKVCFConstants::AS_RAW_RMS_MAPPING_QUALITY_KEY =     "AS_RAW_MQ";
const std::string GATKVCFConstants::AS_CULPRIT_KEY =                     "AS_culprit";
const std::string GATKVCFConstants::AS_VQS_LOD_KEY =                     "AS_VQSLOD";
const std::string GATKVCFConstants::ORIGINAL_AC_KEY =                    "AC_Orig"; //SelectVariants
const std::string GATKVCFConstants::ORIGINAL_AF_KEY =                    "AF_Orig"; //SelectVariants
const std::string GATKVCFConstants::ORIGINAL_AN_KEY =                    "AN_Orig"; //SelectVariants
const std::string GATKVCFConstants::AC_ADJUSTED_KEY =                    "AC_adj"; //GnarlyGenotyper
const std::string GATKVCFConstants::BASE_QUAL_RANK_SUM_KEY =             "BaseQRankSum";
const std::string GATKVCFConstants::BASE_QUAL_HISTOGRAM_KEY =            "BQHIST";
const std::string GATKVCFConstants::AS_BASE_QUAL_RANK_SUM_KEY =          "AS_BaseQRankSum";
const std::string GATKVCFConstants::AS_RAW_BASE_QUAL_RANK_SUM_KEY =      "AS_RAW_BaseQRankSum";
const std::string GATKVCFConstants::GENOTYPE_AND_VALIDATE_STATUS_KEY =   "callStatus";
const std::string GATKVCFConstants::CLIPPING_RANK_SUM_KEY =              "ClippingRankSum";
const std::string GATKVCFConstants::CULPRIT_KEY =                        "culprit";
const std::string GATKVCFConstants::ORIGINAL_DP_KEY =                    "DP_Orig"; //SelectVariants
const std::string GATKVCFConstants::DOWNSAMPLED_KEY =                    "DS";
const std::string GATKVCFConstants::EVENT_COUNT_IN_HAPLOTYPE_KEY =       "ECNT"; //M2
const std::string GATKVCFConstants::FISHER_STRAND_KEY =                  "FS";
const std::string GATKVCFConstants::AS_FISHER_STRAND_KEY =               "AS_FS";
const std::string GATKVCFConstants::AS_SB_TABLE_KEY =                    "AS_SB_TABLE";
const std::string GATKVCFConstants::SB_TABLE_KEY =                       "SB_TABLE";
const std::string GATKVCFConstants::GQ_MEAN_KEY =                        "GQ_MEAN";
const std::string GATKVCFConstants::GQ_STDEV_KEY =                       "GQ_STDDEV";
const std::string GATKVCFConstants::HAPLOTYPE_SCORE_KEY =                "HaplotypeScore";
const std::string GATKVCFConstants::HI_CONF_DENOVO_KEY =                 "hiConfDeNovo";
const std::string GATKVCFConstants::INTERVAL_GC_CONTENT_KEY =            "IGC";
const std::string GATKVCFConstants::INBREEDING_COEFFICIENT_KEY =         "InbreedingCoeff";
const std::string GATKVCFConstants::AS_INBREEDING_COEFFICIENT_KEY =      "AS_InbreedingCoeff";
const std::string GATKVCFConstants::EXCESS_HET_KEY =                     "ExcessHet";
const std::string GATKVCFConstants::RAW_GENOTYPE_COUNT_KEY =             "RAW_GT_COUNT";
const std::string GATKVCFConstants::LIKELIHOOD_RANK_SUM_KEY =            "LikelihoodRankSum";
const std::string GATKVCFConstants::LO_CONF_DENOVO_KEY =                 "loConfDeNovo";
const std::string GATKVCFConstants::MLE_ALLELE_COUNT_KEY =               "MLEAC";
const std::string GATKVCFConstants::MLE_ALLELE_FREQUENCY_KEY =           "MLEAF";
const std::string GATKVCFConstants::MAP_QUAL_RANK_SUM_KEY =              "MQRankSum";
const std::string GATKVCFConstants::RAW_MAP_QUAL_RANK_SUM_KEY =          "RAW_MQRankSum";
const std::string GATKVCFConstants::AS_MAP_QUAL_RANK_SUM_KEY =           "AS_MQRankSum";
const std::string GATKVCFConstants::AS_RAW_MAP_QUAL_RANK_SUM_KEY =       "AS_RAW_MQRankSum";
const std::string GATKVCFConstants::NOCALL_CHROM_KEY =                   "NCC";
const std::string GATKVCFConstants::NUMBER_OF_DISCOVERED_ALLELES_KEY =   "NDA";
const std::string GATKVCFConstants::NEGATIVE_LABEL_KEY =                 "NEGATIVE_TRAIN_SITE";
const std::string GATKVCFConstants::GENOTYPE_PRIOR_KEY =                 "PG";
const std::string GATKVCFConstants::POSITIVE_LABEL_KEY =                 "POSITIVE_TRAIN_SITE";
const std::string GATKVCFConstants::QUAL_BY_DEPTH_KEY =                  "QD";
const std::string GATKVCFConstants::AS_QUAL_BY_DEPTH_KEY =               "AS_QD";
const std::string GATKVCFConstants::AS_QUAL_KEY =                        "AS_QUAL";
const std::string GATKVCFConstants::RAW_QUAL_APPROX_KEY =                "QUALapprox";
const std::string GATKVCFConstants::AS_RAW_QUAL_APPROX_KEY =             "AS_QUALapprox";
const std::string GATKVCFConstants::VARIANT_DEPTH_KEY =                  "VarDP";
const std::string GATKVCFConstants::AS_VARIANT_DEPTH_KEY =               "AS_VarDP";
const std::string GATKVCFConstants::AS_ALT_ALLELE_DEPTH_KEY =            "AS_AltDP";
const std::string GATKVCFConstants::READ_POS_RANK_SUM_KEY =              "ReadPosRankSum";
const std::string GATKVCFConstants::AS_READ_POS_RANK_SUM_KEY =           "AS_ReadPosRankSum";
const std::string GATKVCFConstants::AS_RAW_READ_POS_RANK_SUM_KEY =       "AS_RAW_ReadPosRankSum";
const std::string GATKVCFConstants::REPEATS_PER_ALLELE_KEY =             "RPA";
const std::string GATKVCFConstants::REPEAT_UNIT_KEY =                    "RU";
const std::string GATKVCFConstants::SAMPLE_LIST_KEY =                    "Samples";
const std::string GATKVCFConstants::STRAND_ODDS_RATIO_KEY =              "SOR";
const std::string GATKVCFConstants::AS_STRAND_ODDS_RATIO_KEY =           "AS_SOR";
const std::string GATKVCFConstants::STR_PRESENT_KEY =                    "STR";
const std::string GATKVCFConstants::VQS_LOD_KEY =                        "VQSLOD";
const std::string GATKVCFConstants::CNN_1D_KEY =                         "CNN_1D";
const std::string GATKVCFConstants::CNN_2D_KEY =                         "CNN_2D";
const std::string GATKVCFConstants::F1R2_KEY =                           "F1R2";
const std::string GATKVCFConstants::F2R1_KEY =                           "F2R1";

// Mutect2-specific INFO keys
const std::string GATKVCFConstants::TUMOR_LOG_10_ODDS_KEY =              "TLOD";
const std::string GATKVCFConstants::NORMAL_LOG_10_ODDS_KEY =             "NLOD";
const std::string GATKVCFConstants::IN_PON_KEY =                         "PON";
const std::string GATKVCFConstants::NORMAL_ARTIFACT_LOG_10_ODDS_KEY =    "NALOD";
const std::string GATKVCFConstants::POPULATION_AF_KEY =                  "POPAF";
const std::string GATKVCFConstants::GERMLINE_QUAL_KEY =                  "GERMQ";
const std::string GATKVCFConstants::SEQUENCING_QUAL_KEY =                "SEQQ";
const std::string GATKVCFConstants::POLYMERASE_SLIPPAGE_QUAL_KEY =       "STRQ";
const std::string GATKVCFConstants::STRAND_QUAL_KEY =                    "STRANDQ";
const std::string GATKVCFConstants::CONTAMINATION_QUAL_KEY =             "CONTQ";
const std::string GATKVCFConstants::READ_ORIENTATION_QUAL_KEY =          "ROQ";
const std::string GATKVCFConstants::ORIGINAL_CONTIG_MISMATCH_KEY =       "OCM";
const std::string GATKVCFConstants::N_COUNT_KEY =                        "NCount";
const std::string GATKVCFConstants::AS_UNIQUE_ALT_READ_SET_COUNT_KEY =   "AS_UNIQ_ALT_READ_COUNT";
const std::string GATKVCFConstants::MEDIAN_BASE_QUALITY_KEY =            "MBQ";
const std::string GATKVCFConstants::MEDIAN_MAPPING_QUALITY_KEY =         "MMQ";
const std::string GATKVCFConstants::MEDIAN_FRAGMENT_LENGTH_KEY =         "MFRL";
const std::string GATKVCFConstants::MEDIAN_READ_POSITON_KEY =            "MPOS";
const std::string GATKVCFConstants::UNITIG_SIZES_KEY =                   "UNITIGS";
const std::string GATKVCFConstants::ALIGNMENT_SCORE_DIFFERENCE_KEY =     "ALIGN_DIFF";
const std::string GATKVCFConstants::JOINT_ALIGNMENT_COUNT_KEY =          "NALIGNS";
const std::string GATKVCFConstants::REFERENCE_BASES_KEY =                "REF_BASES";

// Methylation-specific INFO Keys
const std::string GATKVCFConstants::UNCONVERTED_BASE_COVERAGE_KEY =      "UNCONVERTED_BASE_COV";
const std::string GATKVCFConstants::CONVERTED_BASE_COVERAGE_KEY =        "CONVERTED_BASE_COV";
const std::string GATKVCFConstants::METHYLATION_REFERENCE_CONTEXT_KEY =  "REFERENCE_CONTEXT";


// FORMAT keys
const std::string GATKVCFConstants::ALLELE_BALANCE_KEY =                 "AB";
const std::string GATKVCFConstants::JOINT_LIKELIHOOD_TAG_NAME =          "JL"; //FamilyLikelihoodsUtils
const std::string GATKVCFConstants::JOINT_POSTERIOR_TAG_NAME =           "JP"; //FamilyLikelihoodsUtils
const std::string GATKVCFConstants::MIN_DP_FORMAT_KEY =                  "MIN_DP";
const std::string GATKVCFConstants::MAPPING_QUALITY_ZERO_BY_SAMPLE_KEY = "MQ0";
const std::string GATKVCFConstants::HAPLOTYPE_CALLER_PHASING_GT_KEY =    "PGT";
const std::string GATKVCFConstants::HAPLOTYPE_CALLER_PHASING_ID_KEY =    "PID";
const std::string GATKVCFConstants::PHRED_SCALED_POSTERIORS_KEY =        "PP"; //FamilyLikelihoodsUtils / PosteriorLikelihoodsUtils
const std::string GATKVCFConstants::REFERENCE_GENOTYPE_QUALITY =         "RGQ";
const std::string GATKVCFConstants::GENOTYPE_QUALITY_BY_ALLELE_BALANCE = "ABGQ"; //GnarlyGenotyper
const std::string GATKVCFConstants::GENOTYPE_QUALITY_BY_ALT_CONFIDENCE = "ALTGQ"; //GnarlyGenotyper
const std::string GATKVCFConstants::STRAND_COUNT_BY_SAMPLE_KEY =         "SAC";
const std::string GATKVCFConstants::STRAND_BIAS_BY_SAMPLE_KEY =          "SB";
const std::string GATKVCFConstants::FEATURIZED_READ_SETS_KEY =           "FRS";
const std::string GATKVCFConstants::HAPLOTYPE_EQUIVALENCE_COUNTS_KEY =   "HEC";
const std::string GATKVCFConstants::HAPLOTYPE_COMPLEXITY_KEY =           "HAPCOMP";
const std::string GATKVCFConstants::HAPLOTYPE_DOMINANCE_KEY =            "HAPDOM";
const std::string GATKVCFConstants::TRANSMISSION_PROBABILITY_KEY =       "TP"; //PhaseByTransmission
const std::string GATKVCFConstants::FRAGMENT_ALLELE_DEPTHS =             "FAD";

// M2-specific FORMAT keys
const std::string GATKVCFConstants::ALLELE_FRACTION_KEY =                "AF";

//FILTERS
/* Note that many filters used throughout GATK (most notably in VariantRecalibration) are dynamic,
    their names (or descriptions) depend on some threshold.  Those filters are not included here
    */
const std::string GATKVCFConstants::CLUSTERED_EVENTS_FILTER_NAME =                 "clustered_events"; //M2
const std::string GATKVCFConstants::GERMLINE_RISK_FILTER_NAME =                    "germline"; //M2
const std::string GATKVCFConstants::LOW_QUAL_FILTER_NAME =                         "LowQual";
const std::string GATKVCFConstants::ALIGNMENT_ARTIFACT_FILTER_NAME =               "alignment";
const std::string GATKVCFConstants::PON_FILTER_NAME =                              "panel_of_normals"; //M2
const std::string GATKVCFConstants::POLYMERASE_SLIPPAGE =                          "slippage"; //M2
const std::string GATKVCFConstants::TUMOR_EVIDENCE_FILTER_NAME =                   "weak_evidence"; //M2
const std::string GATKVCFConstants::MULTIALLELIC_FILTER_NAME =                     "multiallelic"; //M2
const std::string GATKVCFConstants::STRAND_ARTIFACT_FILTER_NAME =                  "strand_bias"; // M2
const std::string GATKVCFConstants::DUPLICATED_EVIDENCE_FILTER_NAME =              "duplicate";
const std::string GATKVCFConstants::ARTIFACT_IN_NORMAL_FILTER_NAME =               "normal_artifact";
const std::string GATKVCFConstants::MEDIAN_BASE_QUALITY_FILTER_NAME =              "base_qual";
const std::string GATKVCFConstants::MEDIAN_MAPPING_QUALITY_FILTER_NAME =           "map_qual";
const std::string GATKVCFConstants::MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME = "fragment";
const std::string GATKVCFConstants::READ_POSITION_FILTER_NAME =                    "position";
const std::string GATKVCFConstants::CONTAMINATION_FILTER_NAME =                    "contamination";
const std::string GATKVCFConstants::READ_ORIENTATION_ARTIFACT_FILTER_NAME =        "orientation";
const std::string GATKVCFConstants::BAD_HAPLOTYPE_FILTER_NAME =                    "haplotype";
const std::string GATKVCFConstants::STRICT_STRAND_BIAS_FILTER_NAME =               "strict_strand";
const std::string GATKVCFConstants::N_RATIO_FILTER_NAME =                           "n_ratio";
const std::string GATKVCFConstants::ALLELE_FRACTION_FILTER_NAME =                   "low_allele_frac";
const std::string GATKVCFConstants::POSSIBLE_NUMT_FILTER_NAME =                     "possible_numt";
const std::string GATKVCFConstants::LOW_HET_FILTER_NAME =                           "mt_many_low_hets";
const std::string GATKVCFConstants::FAIL =                                           "FAIL";
const std::string GATKVCFConstants::SITE_LEVEL_FILTERS =                             "SITE";


const std::vector<std::string> GATKVCFConstants::MUTECT_FILTER_NAMES = {
    VCFConstants::PASSES_FILTERS_v4, POLYMERASE_SLIPPAGE,
    PON_FILTER_NAME, CLUSTERED_EVENTS_FILTER_NAME, TUMOR_EVIDENCE_FILTER_NAME, GERMLINE_RISK_FILTER_NAME,
    MULTIALLELIC_FILTER_NAME, STRAND_ARTIFACT_FILTER_NAME, ARTIFACT_IN_NORMAL_FILTER_NAME,
    MEDIAN_BASE_QUALITY_FILTER_NAME, MEDIAN_MAPPING_QUALITY_FILTER_NAME,
    MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME,
    READ_POSITION_FILTER_NAME, CONTAMINATION_FILTER_NAME, DUPLICATED_EVIDENCE_FILTER_NAME,
    READ_ORIENTATION_ARTIFACT_FILTER_NAME, BAD_HAPLOTYPE_FILTER_NAME,
    STRICT_STRAND_BIAS_FILTER_NAME, N_RATIO_FILTER_NAME, ALLELE_FRACTION_FILTER_NAME, POSSIBLE_NUMT_FILTER_NAME, FAIL};

const std::vector<std::string> GATKVCFConstants::MUTECT_AS_FILTER_NAMES = {AS_FILTER_STATUS_KEY};

// Symbolic alleles
const std::string GATKVCFConstants::SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG = "ALT";
const std::string GATKVCFConstants::NON_REF_SYMBOLIC_ALLELE_NAME = "NON_REF";
const std::string GATKVCFConstants::ALLELE_SPECIFIC_ANNOTATION_PREFIX = "AS";
const std::string GATKVCFConstants::PSEUDO_DEPTH_KEY = "DD";
const std::string GATKVCFConstants::PSEUDO_FRACTION_KEY = "DF";


