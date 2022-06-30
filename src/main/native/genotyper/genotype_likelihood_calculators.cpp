#include "genotype_likelihood_calculators.hpp"
#include "common/simple_matrix.hpp"
#include "common/utils.hpp"
#include "genotyper/genotype_allele_counts.hpp"
#include "boost/math/special_functions/binomial.hpp"
#include <iostream>
#include <mutex>
#include <vector>



const int GenotypeLikelihoodCalculators::MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY = 500500;

const int GenotypeLikelihoodCalculators::GENOTYPE_COUNT_OVERFLOW = -1;


SimpleMatrix<int> 
GenotypeLikelihoodCalculators::buildAlleleFirstGenotypeOffsetTable(
    int max_ploidy, int max_allele)
{
    int row_count = max_ploidy + 1;
    int col_count = max_allele + 1;

    SimpleMatrix<int> result(row_count, col_count);

    result.fill(0, result.rows(), 0, 1, 0);
    result.fill(0, result.rows(), 1, result.cols(), 1);

    for (int ploidy = 1; ploidy < row_count; ploidy++) {
        for (int allele = 1; allele < col_count; allele++) {
            result(ploidy, allele) =
                result(ploidy, allele - 1) + result(ploidy - 1, allele);
            if (result(ploidy, allele) > MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY  // it is weired that GATK not use this line to check genotype count overflow
                || result(ploidy, allele) < result(ploidy, allele - 1))
            {
                result(ploidy, allele) = GENOTYPE_COUNT_OVERFLOW;
            }
        }
    }
    return result;
}

std::vector<std::vector<GenotypeAlleleCounts>>
GenotypeLikelihoodCalculators::buildGenotypeAlleleCountsTable(
    int max_ploidy, int max_allele,
    const SimpleMatrix<int> &genotypeOffsetTable)
{
    const int row_counts = max_ploidy + 1;
    std::vector<std::vector<GenotypeAlleleCounts>> result(row_counts);
    for (int ploidy = 0; ploidy < row_counts; ploidy++) {
        result[ploidy] = buildGenotypeAlleleCountsArray(
            ploidy, max_allele, genotypeOffsetTable);
    }

    return result;
}

std::vector<GenotypeAlleleCounts>
GenotypeLikelihoodCalculators::buildGenotypeAlleleCountsArray(
    int ploidy, int allele_count, const SimpleMatrix<int> &genotypeOffsetTable)
{
    const int length = genotypeOffsetTable(ploidy, allele_count);
    const int strong_ref_length = length == GENOTYPE_COUNT_OVERFLOW ?
        MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY :
        std::min(length, MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY);
    std::vector<GenotypeAlleleCounts> result(strong_ref_length);
    result[0] = GenotypeAlleleCounts::first(ploidy);
    for (int i = 1; i < strong_ref_length; ++i) {
        result[i] = result[i-1].next();
    }
    return result;    
}


void GenotypeLikelihoodCalculators::ensureCapacity(
    int requested_max_ploidy, int requested_max_allele)
{
    const bool needs_to_expand_allele_capacity =
        requested_max_allele > max_allele_;
    const bool needs_to_expand_ploidy_capacity =
        requested_max_ploidy > max_ploidy_;
    
    // no need to expand tables
    if (!needs_to_expand_allele_capacity && !needs_to_expand_ploidy_capacity) {
        return;
    }

    max_ploidy_ = std::max(max_ploidy_, requested_max_ploidy);
    max_allele_ = std::max(max_allele_, requested_max_allele);

    // update allele_first_genotype_offset_by_ploidy_ table
    allele_first_genotype_offset_by_ploidy_ =
        buildAlleleFirstGenotypeOffsetTable(max_ploidy_, max_allele_);
    
    genotype_table_by_ploidy_ = buildGenotypeAlleleCountsTable(
        max_ploidy_, max_allele_, allele_first_genotype_offset_by_ploidy_);
}


int GenotypeLikelihoodCalculators::genotypeCount(
    int ploidy, int allele_count) const
{
    const int result = allele_first_genotype_offset_by_ploidy_(
        ploidy, allele_count);
    if (result == GENOTYPE_COUNT_OVERFLOW) {
        double large_genotype_count = boost::math::binomial_coefficient<double>(
            ploidy + allele_count - 1, allele_count - 1);
        std::cerr << "the number of genotypes is too large for ploidy:"
            << ploidy << " and allele count:" << allele_count
            << ". number of genotypes:" << large_genotype_count;
        std::exit(1);
    }
    return result;
}


void
GenotypeLikelihoodCalculators::calculateGenotypeCountsUsingTablesAndValidate(
    int ploidy, int allele_count)
{
    if (calculateGenotypeCountUsingTables(ploidy, allele_count) ==
        GENOTYPE_COUNT_OVERFLOW)
    {
        double large_genotype_count = boost::math::binomial_coefficient<double>(
            ploidy + allele_count - 1, allele_count - 1);
        std::cerr << "the number of genotypes is too large for ploidy:"
            << ploidy << " and allele count:" << allele_count
            << ". number of genotypes:" << large_genotype_count;
        std::exit(1);
    }
}


void GenotypeLikelihoodCalculators::checkPloidyAndMaximumAllele(
    int ploidy, int max_allele)
{
    Utils::validateArg(ploidy, "the ploidy can not be negative");
    Utils::validateArg(max_allele, "the max_allele can not be negative");
}

void GenotypeLikelihoodCalculators::checkOffsetTableCapacity(
    const SimpleMatrix<int> &offset_table, int max_ploidy, int max_allele)
{
    Utils::validateArg(offset_table.rows() > max_ploidy, "the allele first "
        "genotype offset by ploidy table not have enough rows for ploidy.");
    Utils::validateArg(offset_table.cols() > max_allele, "the allele first "
        "genotype offset by ploidy table not have enough cols for allele.");
}


int GenotypeLikelihoodCalculators::calculateGenotypeCountUsingTables(
    int ploidy, int allele_count)
{
    checkPloidyAndMaximumAllele(ploidy, allele_count);
    ensureCapacity(ploidy, allele_count);
    return allele_first_genotype_offset_by_ploidy_(ploidy, allele_count);
}
