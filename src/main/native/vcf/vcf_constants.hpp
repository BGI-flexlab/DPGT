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
#ifndef DPGT_VCF_CONSTANTS_HPP
#define DPGT_VCF_CONSTANTS_HPP

#include <string>
#include <sys/cdefs.h>


class VCFConstants {
public:
    // static const Locale VCF_LOCALE;

    // reserved INFO/FORMAT field keys
    static const std::string ANCESTRAL_ALLELE_KEY;
    static const std::string ALLELE_COUNT_KEY;
    static const std::string ALLELE_FREQUENCY_KEY;
    static const std::string ALLELE_NUMBER_KEY;
    static const std::string RMS_BASE_QUALITY_KEY;
    static const std::string CIGAR_KEY;
    static const std::string DBSNP_KEY;
    static const std::string DEPTH_KEY;
    static const std::string END_KEY;

    static const std::string GENOTYPE_FILTER_KEY;
    static const std::string GENOTYPE_KEY;
    static const std::string GENOTYPE_POSTERIORS_KEY;
    static const std::string GENOTYPE_QUALITY_KEY;
    static const std::string GENOTYPE_ALLELE_DEPTHS;
    static const std::string GENOTYPE_PL_KEY;
    static const std::string EXPECTED_ALLELE_COUNT_KEY;
    // deprecated
    static const std::string GENOTYPE_LIKELIHOODS_KEY;

    static const std::string HAPMAP2_KEY;
    static const std::string HAPMAP3_KEY;
    static const std::string HAPLOTYPE_QUALITY_KEY;
    static const std::string RMS_MAPPING_QUALITY_KEY;
    static const std::string MAPPING_QUALITY_ZERO_KEY;
    static const std::string SAMPLE_NUMBER_KEY;
    static const std::string PHASE_QUALITY_KEY;
    static const std::string PHASE_SET_KEY;
    static const std::string OLD_DEPTH_KEY;
    static const std::string STRAND_BIAS_KEY;
    static const std::string SOMATIC_KEY;
    static const std::string VALIDATED_KEY;
    static const std::string THOUSAND_GENOMES_KEY;
    
    // reserved INFO for structural variants
    /** INFO Type of structural variant */
    static const std::string SVTYPE;

    // separators
    static const std::string FORMAT_FIELD_SEPARATOR;
    static const std::string GENOTYPE_FIELD_SEPARATOR;
    static const char   GENOTYPE_FIELD_SEPARATOR_CHAR;
    static const std::string FIELD_SEPARATOR;
    static const char   FIELD_SEPARATOR_CHAR;
    static const std::string FILTER_CODE_SEPARATOR;
    static const std::string INFO_FIELD_ARRAY_SEPARATOR;
    static const char INFO_FIELD_ARRAY_SEPARATOR_CHAR;
    static const std::string ID_FIELD_SEPARATOR;
    static const std::string INFO_FIELD_SEPARATOR;
    static const char INFO_FIELD_SEPARATOR_CHAR;
    static const std::string UNPHASED;
    static const std::string PHASED;
    static const std::string PHASED_SWITCH_PROB_v3;
    static const std::string PHASING_TOKENS;

    // header lines
    static const std::string FILTER_HEADER_START;
    static const std::string FORMAT_HEADER_START;
    static const std::string INFO_HEADER_START;
    static const std::string ALT_HEADER_KEY;
    static const std::string ALT_HEADER_START;
    static const std::string CONTIG_HEADER_KEY;
    static const std::string CONTIG_HEADER_START;

    static const int ALT_HEADER_OFFSET;

    static const std::string PEDIGREE_HEADER_KEY;
    static const std::string PEDIGREE_HEADER_START;
    static const int PEDIGREE_HEADER_OFFSET;

    static const std::string SAMPLE_HEADER_KEY;
    static const std::string SAMPLE_HEADER_START;
    static const int SAMPLE_HEADER_OFFSET;

    static const std::string META_HEADER_KEY;
    static const std::string META_HEADER_START;
    static const int META_HEADER_OFFSET;

    // old indel alleles
    static const char DELETION_ALLELE_v3;
    static const char INSERTION_ALLELE_v3;

    // special alleles
    static const char SPANNING_DELETION_ALLELE;
    static const char NO_CALL_ALLELE;
    static const char NULL_ALLELE;


    // missing/default values
    static const std::string UNFILTERED;
    static const std::string PASSES_FILTERS_v3;
    static const std::string PASSES_FILTERS_v4;
    static const std::string EMPTY_ID_FIELD;
    static const std::string EMPTY_INFO_FIELD;
    static const std::string EMPTY_ALTERNATE_ALLELE_FIELD;
    static const std::string MISSING_VALUE_v4;
    static const std::string MISSING_QUALITY_v3;
    static const double MISSING_QUALITY_v3_DOUBLE;

    static const std::string MISSING_GENOTYPE_QUALITY_v3;
    static const std::string MISSING_HAPLOTYPE_QUALITY_v3;
    static const std::string MISSING_DEPTH_v3;
    static const std::string UNBOUNDED_ENCODING_v4;
    static const std::string UNBOUNDED_ENCODING_v3;
    static const std::string PER_ALTERNATE_COUNT;
    static const std::string PER_ALLELE_COUNT;
    static const std::string PER_GENOTYPE_COUNT;
    static const std::string EMPTY_ALLELE;
    static const std::string EMPTY_GENOTYPE;
    static const int MAX_GENOTYPE_QUAL;

    static const double VCF_ENCODING_EPSILON;
};

#endif  // DPGT_VCF_CONSTANTS_HPP
