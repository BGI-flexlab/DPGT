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
#include "genotype_likelihood_calculator.hpp"
#include "common/utils.hpp"
#include "genotyper/genotype_allele_counts.hpp"
#include "genotyper/genotype_likelihood_calculators.hpp"
#include <cstdlib>
#include <functional>
#include <iostream>
#include <queue>
#include <tuple>
#include <vector>


GenotypeLikelihoodCalculator::GenotypeLikelihoodCalculator(
    int ploidy, int allele_count,
    const SimpleMatrix<int> *allele_first_genotype_offset_by_ploidy,
    const std::vector<GenotypeAlleleCounts> *genotype_allele_counts):
    ploidy_(ploidy), allele_count_(allele_count),
    allele_first_genotype_offset_by_ploidy_(allele_first_genotype_offset_by_ploidy),
    genotype_allele_counts_(genotype_allele_counts)
{
    max_distinct_alleles_in_genotype_ = std::max(ploidy_, allele_count_);
    genotype_count_ = (*allele_first_genotype_offset_by_ploidy)(
        ploidy_, allele_count_);
}

int GenotypeLikelihoodCalculator::allelesToIndex(
    const std::vector<int> &allele_count_array)
{
    Utils::validateArg(allele_count_array.size() > 0,
        "allele count array can not be empty");
    Utils::validateArg((allele_count_array.size() & 1) == 0,
        "allele count array can not have odd size");
    // clear the allele heap
    allele_heap_ =
        std::priority_queue<int, std::vector<int>>();
    int index, count;
    for (int i = 0; i < static_cast<int>(allele_count_array.size()); i+=2) {
        index = allele_count_array[i];
        count = allele_count_array[i+1];
        Utils::validateArg(count >= 0, "allele can not be less than 0");
        for (int j = 0; j < count; ++j) {
            allele_heap_.push(index);
        }
    }
    return alleleHeapToIndex();
}


int GenotypeLikelihoodCalculator::alleleHeapToIndex() {
    Utils::validateArg(static_cast<int>(allele_heap_.size()) == ploidy_,
        "the sum of allele counts must be equal to the ploidy of the calculator");
    Utils::validateArg(allele_heap_.top() < allele_count_,
        "allele index must less than allele count");
    int result = 0;
    int allele;
    for (int p = ploidy_; p > 0; --p) {
        allele = allele_heap_.top();
        allele_heap_.pop();
        result += (*allele_first_genotype_offset_by_ploidy_)(p, allele);
    }
    return result;
}


int *GenotypeLikelihoodCalculator::genotypeIndexMap(
    int *old_to_new_allele_index_map, int new_allele_count,
    const GenotypeLikelihoodCalculators &calculators, int &n)
{
    Utils::validateArg(new_allele_count <= allele_count_, "this calculator "
        "dose not have enough capacity for input allele count");
    n = new_allele_count == allele_count_ ?
        genotype_count_ :
        calculators.genotypeCount(ploidy_, new_allele_count);
    
    int *result = (int *)calloc(n, sizeof(int));
    
    // buffer 
    std::vector<int> sorted_allele_counts(std::max(ploidy_, allele_count_) << 1);
    // clear the allele heap
    allele_heap_ =
        std::priority_queue<int, std::vector<int>>();
    GenotypeAlleleCounts allele_counts = (*genotype_allele_counts_)[0];
    for (int i = 0; i < n; ++i) {
        genotypeIndexMapPerGenotypeIndex(i, allele_counts,
            old_to_new_allele_index_map, result, sorted_allele_counts);
        allele_counts = nextGenotypeAlleleCounts(allele_counts);
    }
    return result;
}


GenotypeAlleleCounts GenotypeLikelihoodCalculator::nextGenotypeAlleleCounts(
        const GenotypeAlleleCounts &allele_counts) const
{
    const int index = allele_counts.index(); // genotype index
    GenotypeAlleleCounts result;
    const int next = index + 1;
    if (next < static_cast<int>(genotype_allele_counts_->size())) {
        result = genotype_allele_counts_->at(next);
    } else  {  // next genotype index out of cache range, calculate next
        result = allele_counts;
        result.increase();
    }
    return result;
}


void GenotypeLikelihoodCalculator::genotypeIndexMapPerGenotypeIndex(
    int new_genotype_index,
    const GenotypeAlleleCounts &allele_counts,
    int *old_to_new_allele_index_map,
    int *destination,
    std::vector<int> &sorted_allele_counts_buffer)
{
    const int distinct_allele_count = allele_counts.distinctAlleleCount();
    allele_counts.copyAlleleCounts(sorted_allele_counts_buffer, 0);
    for (int j = 0, jj = 0; j < distinct_allele_count; ++j) {
        const int new_index = sorted_allele_counts_buffer[jj++];
        const int repeats = sorted_allele_counts_buffer[jj++];
        const int old_index = old_to_new_allele_index_map[new_index];
        if (old_index < 0 || new_index > allele_count_) {
            std::cerr << "find invalid old allele index: " << old_index << " for"
                " new allele index: " << new_index << std::endl;
            std::exit(1);
        }
        for (int k = 0; k < repeats; ++k) {
            allele_heap_.push(old_index);
        }
    }
    destination[new_genotype_index] = alleleHeapToIndex();
}

GenotypeAlleleCounts
GenotypeLikelihoodCalculator::genotypeAlleleCountsAt(int index)
{
    if (index < static_cast<int>(genotype_allele_counts_->size())) {
        return (*genotype_allele_counts_)[index];
    } else {
        // request index is outof lookup table capacity, calculate it
        GenotypeAlleleCounts result = genotype_allele_counts_->back();
        result.increase(index - genotype_allele_counts_->size() + 1);
        return result;
    }
}


GenotypeLikelihoodCalculator createGenotypeLikelihoodCalculator(int ploidy,
    int allele_count, GenotypeLikelihoodCalculators &calculators)
{
    calculators.calculateGenotypeCountsUsingTablesAndValidate(
        ploidy, allele_count);
    return GenotypeLikelihoodCalculator(ploidy, allele_count,
        &(calculators.alleleFirstGenotypeOffsetByPloidy()),
        &(calculators.genotypeTableByPloidy()[ploidy]));
}
