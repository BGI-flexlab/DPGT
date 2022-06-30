#ifndef DPGT_GENOTYPE_LIKELIHOOD_CALCULATORS_HPP
#define DPGT_GENOTYPE_LIKELIHOOD_CALCULATORS_HPP

#include "common/simple_matrix.hpp"
#include "genotyper/genotype_allele_counts.hpp"
#include <vector>
#include <mutex>


/**
 * @brief Genotype likelihood calculator utility(lookup tables).
 * 
 * This class provide genotype likelihood calculators with any number of alleles
 * able given an arbitrary ploidy and allele count (number of distinct alleles).
 * 
 */
class GenotypeLikelihoodCalculators {

public:
    // maximum number of genotypes that this caclulator can handle
    static const int  MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY;

    // Mark to indicate genotype-count overflow due to a large number of
    // allele and ploidy
    static const int GENOTYPE_COUNT_OVERFLOW;

private:
    int max_ploidy_ = 2;  // maximum ploidy supported by the tables
    int max_allele_ = 1;  // maximum allele supported by the tables

    /**
     * @brief table of first genotype offset given ploidy and allele index
     * row is ploidy, column is allele index
     */
    SimpleMatrix<int> allele_first_genotype_offset_by_ploidy_;

    /**
     * @brief table of first genotype(GenotypeAlleleCounts) given ploidy
     * and allele index
     * row is ploidy, column is allele index
     */
    std::vector<std::vector<GenotypeAlleleCounts>> genotype_table_by_ploidy_;

public:
    GenotypeLikelihoodCalculators() {
        allele_first_genotype_offset_by_ploidy_ =
            buildAlleleFirstGenotypeOffsetTable(max_ploidy_, max_allele_);
        buildGenotypeAlleleCountsTable(max_ploidy_, max_allele_,
            allele_first_genotype_offset_by_ploidy_);
    }



    /**
     * @brief build table of first genotype offset by ploidy and allele index
     *     
     *     This value are calculated recursively as follows:
     *
     *     Offset[p][a] := Offset[p-1][a] + Offset[p][a-1] when a > 0, p > 0
     *                     0                               when a == 0
     *                     1                               otherwise
     *
     *
     *         0 1 1  1  1  1   1 ...
     *         0 1 2  3  4  5   6 ...
     *         0 1 3  6 10 15  21 ...
     *         0 1 4 10 20 35  56 ...
     *         0 1 5 15 35 70 126 ...
     *         0 ..................
     *    
     *    Offset[p][a] is equivalent to number of genotypes given ploidy and
     *    number of alleles
     *    Offset[p][a] == ${P+A-1 \choose P}$
     * 
     * @param max_ploidy_ maximum ploidy for the table
     * @param max_allele_ maximum allele for the table
     * @return SimpleMatrix<int> table of first genotype offset
     * by ploidy and allele index, result have shape (max_ploidy_+1, max_allele_+1)
     */
    static SimpleMatrix<int> buildAlleleFirstGenotypeOffsetTable(
        int max_ploidy, int max_allele);
    


    static std::vector<std::vector<GenotypeAlleleCounts>>
    buildGenotypeAlleleCountsTable(int max_ploidy, int max_allele,
        const SimpleMatrix<int> &genotypeOffsetTable);

    /**
     * @brief build a genotype allele counts array given the genotype ploidy and
     * how many genotypes you need
     * 
     * @param ploidy requested ploidy
     * @param allele_count number of different alleles that the genotype table
     * must support
     * @param genotypeOffsetTable table of first genotype offset
     * by ploidy and allele index
     * @return std::vector<GenotypeAlleleCounts> 
     */
    static std::vector<GenotypeAlleleCounts> buildGenotypeAlleleCountsArray(
        int ploidy, int allele_count,
        const SimpleMatrix<int> &genotypeOffsetTable);
    
    /**
     * @brief get genotype count given ploidy and allele count
     * 
     * @param ploidy input ploidy
     * @param allele_count input allele count
     * @return int genotype count
     */
    int genotypeCount(int ploidy, int allele_count) const;

    /**
     * @brief Calculate genotype counts using the tables and validate that
     * there is no overflow
     * call this function will expand lookup tables if necessary
     * 
     * @param ploidy input ploidy
     * @param allele_count input allele count
     */
    void calculateGenotypeCountsUsingTablesAndValidate(int ploidy,
        int allele_count);
    
    const SimpleMatrix<int> &alleleFirstGenotypeOffsetByPloidy() const {
        return allele_first_genotype_offset_by_ploidy_;
    }

    const std::vector<std::vector<GenotypeAlleleCounts>> &
    genotypeTableByPloidy() const {
        return genotype_table_by_ploidy_;
    }

private:
    /**
     * @brief expand tables as need
     * 
     * @param requested_max_ploidy the new requested max ploidy
     * @param requested_max_allele the new requested max allele
     */
    void ensureCapacity(int requested_max_ploidy, int requested_max_allele);


    static void checkPloidyAndMaximumAllele(int ploidy, int max_allele);

    static void checkOffsetTableCapacity(
        const SimpleMatrix<int> &offset_table, int max_ploid, int max_allele);

    int calculateGenotypeCountUsingTables(int ploidy, int allele_count);
};



#endif  // DPGT_GENOTYPE_LIKELIHOOD_CALCULATORS_HPP
