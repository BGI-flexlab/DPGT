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
#ifndef DPGT_GENOTYPE_LIKELIHOOD_CACULATOR_HPP
#define DPGT_GENOTYPE_LIKELIHOOD_CACULATOR_HPP


#include "common/simple_matrix.hpp"
#include "genotyper/genotype_allele_counts.hpp"
#include "genotyper/genotype_likelihood_calculators.hpp"
#include <functional>
#include <vector>
#include <queue>

class GenotypeLikelihoodCalculator {
private:
    /**
     * @brief ploidy for this calculator
     */
    int ploidy_;

    /**
     * @brief number of alleles for this calculator
     */
    int allele_count_;

    /**
     * @brief offset table for this calculator
     * pointer to GenotypeLikelihoodCalculator allele_first_genotype_offset_by_ploidy_
     * 
     * allele_first_genotype_offset_by_ploidy_ in GenotypeLikelihoodCalculator
     * is used as a cached table
     * 
     */
    const SimpleMatrix<int> *allele_first_genotype_offset_by_ploidy_ = nullptr;

    /**
     * @brief genotype table for this calculator 
     */
    const std::vector<GenotypeAlleleCounts> *genotype_allele_counts_ = nullptr;

    /**
     * @brief number of genotype for this calculator(ploidy and allele count)
     */
    int genotype_count_;

    /**
     * @brief maximum distinct alleles for this calculator
     */
    int max_distinct_alleles_in_genotype_;

    /**
     * @brief max-heap used for calculate genotype index from allele indices
     */
    std::priority_queue<int, std::vector<int>> allele_heap_;

    // TODO caculate genotype likelihoods

public:

    /**
     * @brief Construct a new Genotype Likelihood Calculator object
     * 
     * @param ploidy 
     * @param allele_count number of alleles
     * @param allele_first_genotype_offset_by_ploidy_ offset table for this calculator
     * @param genotype_allele_counts_ genotype table for this calculator
     */
    GenotypeLikelihoodCalculator(int ploidy, int allele_count,
        const SimpleMatrix<int> *allele_first_genotype_offset_by_ploidy,
        const std::vector<GenotypeAlleleCounts> *genotype_allele_counts);
    
    /**
     * @brief calculate genotype index given allele indices of a genotype
     * 
     * @param allele_count_array genotype represent by an allele count array as 
     * described in genotype_allele_counts.cpp
     * @return int genotype index
     */
    int allelesToIndex(const std::vector<int> &allele_count_array);

    /**
     * @brief calculate new genotype index map to old genotype index
     * given old to new allele index map
     * 
     * @param old_to_new_allele_index_map old to new allele index map, vector 
     * index is new allele index, vector value is corresponding old allele index
     * @param calculators GenotypeLikelihoodCalculators
     * @param n result length
     * @return result genotype index map, offset is new genotype index, value is
     * old genotype index
     */
    int *genotypeIndexMap(
        int *old_to_new_allele_index_map, int new_allele_count,
        const GenotypeLikelihoodCalculators &calculators,
        int &n);

    /**
     * @return int ploidy of this calculator
     */
    int ploidy() const {
        return ploidy_;
    }

    /**
     * @return int allele count of this calculator
     */
    int alleleCount() const {
        return allele_count_;
    }

    /**
     * @return int genotype count of this calculator
     */
    int genotypeCount() const {
        return genotype_count_;
    }
    
    /**
     * @return int maximum distinct alleles for this calculator
     */
    int maxDistinctAlleleInGenotype() const {
        return max_distinct_alleles_in_genotype_;
    }

    /**
     * @brief get genotype allele counts(represent a genotype) by PL index
     * 
     * @param index PL index
     * @return GenotypeAlleleCounts genotype represent by genotype allele counts
     */
    GenotypeAlleleCounts genotypeAlleleCountsAt(int index);

private:

    /**
     * @brief calculate genotype index from the internal allele-heap
     * 
     * @return int genotype index
     */
    int alleleHeapToIndex();

    GenotypeAlleleCounts nextGenotypeAlleleCounts(
        const GenotypeAlleleCounts &allele_counts) const;

    void genotypeIndexMapPerGenotypeIndex(int new_genotype_index,
        const GenotypeAlleleCounts &allele_counts,
        int *old_to_new_allele_index_map,
        int *destination,
        std::vector<int> &sorted_allele_counts_buffer);
};

/**
 * @brief Create a Genotype Likelihood Calculator object
 *
 * call this function will change the calculators if ploidy larger
 * than maximum ploidy of the calculators or allele count larger maximum allele
 * count of the calculators
 * 
 * @param ploidy ploidy
 * @param allele_count number of alleles
 * @param calculators Genotype Likelihood Calculator lookup tables
 * @return GenotypeLikelihoodCalculator
 */
GenotypeLikelihoodCalculator createGenotypeLikelihoodCalculator(int ploidy,
    int allele_count, GenotypeLikelihoodCalculators &calculators);



#endif  // DPGT_GENOTYPE_LIKELIHOOD_CACULATOR_HPP
