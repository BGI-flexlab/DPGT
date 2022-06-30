#include "vcf_constants.hpp"


// const Locale VCF_LOCALE = Locale.US;

// reserved INFO/FORMAT field keys
const std::string VCFConstants::ANCESTRAL_ALLELE_KEY = "AA";
const std::string VCFConstants::ALLELE_COUNT_KEY = "AC";
const std::string VCFConstants::ALLELE_FREQUENCY_KEY = "AF";
const std::string VCFConstants::ALLELE_NUMBER_KEY = "AN";
const std::string VCFConstants::RMS_BASE_QUALITY_KEY = "BQ";
const std::string VCFConstants::CIGAR_KEY = "CIGAR";
const std::string VCFConstants::DBSNP_KEY = "DB";
const std::string VCFConstants::DEPTH_KEY = "DP";
const std::string VCFConstants::END_KEY = "END";

const std::string VCFConstants::GENOTYPE_FILTER_KEY = "FT";
const std::string VCFConstants::GENOTYPE_KEY = "GT";
const std::string VCFConstants::GENOTYPE_POSTERIORS_KEY = "GP";
const std::string VCFConstants::GENOTYPE_QUALITY_KEY = "GQ";
const std::string VCFConstants::GENOTYPE_ALLELE_DEPTHS = "AD"; //AD isn't reserved, but is specifically handled by VariantContext
const std::string VCFConstants::GENOTYPE_PL_KEY = "PL";   // phred-scaled genotype likelihoods
const std::string VCFConstants::EXPECTED_ALLELE_COUNT_KEY = "EC";

const std::string VCFConstants::GENOTYPE_LIKELIHOODS_KEY = "GL";         // log10 scaled genotype likelihoods

const std::string VCFConstants::HAPMAP2_KEY = "H2";
const std::string VCFConstants::HAPMAP3_KEY = "H3";
const std::string VCFConstants::HAPLOTYPE_QUALITY_KEY = "HQ";
const std::string VCFConstants::RMS_MAPPING_QUALITY_KEY = "MQ";
const std::string VCFConstants::MAPPING_QUALITY_ZERO_KEY = "MQ0";
const std::string VCFConstants::SAMPLE_NUMBER_KEY = "NS";
const std::string VCFConstants::PHASE_QUALITY_KEY = "PQ";
const std::string VCFConstants::PHASE_SET_KEY = "PS";
const std::string VCFConstants::OLD_DEPTH_KEY = "RD";
const std::string VCFConstants::STRAND_BIAS_KEY = "SB";
const std::string VCFConstants::SOMATIC_KEY = "SOMATIC";
const std::string VCFConstants::VALIDATED_KEY = "VALIDATED";
const std::string VCFConstants::THOUSAND_GENOMES_KEY = "1000G";

// reserved INFO for structural variants
/** INFO Type of structural variant */
const std::string VCFConstants::SVTYPE = "SVTYPE";    

// separators
const std::string VCFConstants::FORMAT_FIELD_SEPARATOR = ":";
const std::string VCFConstants::GENOTYPE_FIELD_SEPARATOR = ":";
const char VCFConstants::GENOTYPE_FIELD_SEPARATOR_CHAR = ':';
const std::string VCFConstants::FIELD_SEPARATOR = "\t";
const char VCFConstants::FIELD_SEPARATOR_CHAR = '\t';
const std::string VCFConstants::FILTER_CODE_SEPARATOR = ";";
const std::string VCFConstants::INFO_FIELD_ARRAY_SEPARATOR = ",";
const char VCFConstants::INFO_FIELD_ARRAY_SEPARATOR_CHAR= ',';
const std::string VCFConstants::ID_FIELD_SEPARATOR = ";";
const std::string VCFConstants::INFO_FIELD_SEPARATOR = ";";
const char VCFConstants::INFO_FIELD_SEPARATOR_CHAR = ';';
const std::string VCFConstants::UNPHASED = "/";
const std::string VCFConstants::PHASED = "|";
const std::string VCFConstants::PHASED_SWITCH_PROB_v3 = "\\";
const std::string VCFConstants::PHASING_TOKENS = "/|\\";

// header lines
const std::string VCFConstants::FILTER_HEADER_START = "##FILTER";
const std::string VCFConstants::FORMAT_HEADER_START = "##FORMAT";
const std::string VCFConstants::INFO_HEADER_START = "##INFO";
const std::string VCFConstants::ALT_HEADER_KEY = "ALT";
const std::string VCFConstants::ALT_HEADER_START = "##" + ALT_HEADER_KEY ;
const std::string VCFConstants::CONTIG_HEADER_KEY = "contig";
const std::string VCFConstants::CONTIG_HEADER_START = "##" + CONTIG_HEADER_KEY;

const int VCFConstants::ALT_HEADER_OFFSET = ALT_HEADER_START.length() + 1;

const std::string VCFConstants::PEDIGREE_HEADER_KEY = "PEDIGREE";
const std::string VCFConstants::PEDIGREE_HEADER_START = "##" + PEDIGREE_HEADER_KEY;
const int VCFConstants::PEDIGREE_HEADER_OFFSET = PEDIGREE_HEADER_START.length() + 1;

const std::string VCFConstants::SAMPLE_HEADER_KEY = "SAMPLE";
const std::string VCFConstants::SAMPLE_HEADER_START = "##" + SAMPLE_HEADER_KEY;
const int VCFConstants::SAMPLE_HEADER_OFFSET = SAMPLE_HEADER_START.length() + 1;

const std::string VCFConstants::META_HEADER_KEY = "META";
const std::string VCFConstants::META_HEADER_START = "##" + META_HEADER_KEY;
const int VCFConstants::META_HEADER_OFFSET = META_HEADER_START.length() + 1;

// old indel alleles
const char VCFConstants::DELETION_ALLELE_v3 = 'D';
const char VCFConstants::INSERTION_ALLELE_v3 = 'I';

// special alleles
const char VCFConstants::SPANNING_DELETION_ALLELE = '*';
const char VCFConstants::NO_CALL_ALLELE = '.';
const char VCFConstants::NULL_ALLELE = '-';


// missing/default values
const std::string VCFConstants::UNFILTERED = ".";
const std::string VCFConstants::PASSES_FILTERS_v3 = "0";
const std::string VCFConstants::PASSES_FILTERS_v4 = "PASS";
const std::string VCFConstants::EMPTY_ID_FIELD = ".";
const std::string VCFConstants::EMPTY_INFO_FIELD = ".";
const std::string VCFConstants::EMPTY_ALTERNATE_ALLELE_FIELD = ".";
const std::string VCFConstants::MISSING_VALUE_v4 = ".";
const std::string VCFConstants::MISSING_QUALITY_v3 = "-1";
const double VCFConstants::MISSING_QUALITY_v3_DOUBLE = -1.0;

const std::string VCFConstants::MISSING_GENOTYPE_QUALITY_v3 = "-1";
const std::string VCFConstants::MISSING_HAPLOTYPE_QUALITY_v3 = "-1";
const std::string VCFConstants::MISSING_DEPTH_v3 = "-1";
const std::string VCFConstants::UNBOUNDED_ENCODING_v4 = ".";
const std::string VCFConstants::UNBOUNDED_ENCODING_v3 = "-1";
const std::string VCFConstants::PER_ALTERNATE_COUNT = "A";
const std::string VCFConstants::PER_ALLELE_COUNT = "R";
const std::string VCFConstants::PER_GENOTYPE_COUNT = "G";
const std::string VCFConstants::EMPTY_ALLELE = ".";
const std::string VCFConstants::EMPTY_GENOTYPE = "./.";
const int VCFConstants::MAX_GENOTYPE_QUAL = 99;

const double VCFConstants::VCF_ENCODING_EPSILON = 0.00005; // when we consider fields equal(), used in the Qual compare


