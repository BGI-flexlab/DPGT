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
#ifndef DPGT_GENOTYPE_ALLELE_COUNTS_HPP
#define DPGT_GENOTYPE_ALLELE_COUNTS_HPP

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "boost/algorithm/string/join.hpp"
#include "boost/math/constants/constants.hpp"
#include "common/math_utils.hpp"
#include "vcf/vcf_constants.hpp"


/**
 * Collection of allele counts for a genotype. It encompasses what alleles are
 * present in the genotype and in what number.
 *
 * Genotypes are represented as a single array of alternating alleles and counts
 * where only alleles with non-zero counts are included:
 * [allele1, count1, allele2, count2...]
 * 
 * allele is represent by it's index(rank) in all alleles of the variant.
 * 
 * For example:
 * [0, 1, 2, 1] has two alleles with indices 0 and 2, both with count 1,
 * corresponding to diploid genotype 0/2.
 * 
 * [2, 1, 4, 2, 7, 1] has three alleles with indices 2, 4 and 7 have count 1, 2
 * and 1 respectively. It is corresponds to tetraploid genotype 2/4/4/7
 * 
 * [0, 0, 1, 2] is not valid because allele 0 has count of 0 and should be absent
 * from the array.
 * [1, 1, 0, 1] is not valid because allele 1 comes before allele 0.
 * 
 * instances have themselves their own index (returned by {@link #index() index()},
 * that indicate their 0-based ordinal within the possible genotype combinations
 * with the same ploidy
 *
 * The total number of possible genotypes is only bounded by the maximum allele
 * index.
 */
class GenotypeAlleleCounts {
private:
    static const double UNCOMPUTED_LOG_10_COMBINATION_COUNT;

    /**
     * The log10 number of phased genotypes corresponding to this unphased genotype.
     * For example,
     * [0, 1, 1, 1] = AB:  log10(2)
     * [0, 2] = AA:  log10(1)
     * [0, 1, 1, 1, 2, 1] = ABC: log10(6)
     * [0, 2, 1, 2] = AABB: log10(4!/(2!2!))
     * This is evaluated lazily and only calculated if its getter is invoked.
     */
    double log10_combination_count_ = UNCOMPUTED_LOG_10_COMBINATION_COUNT;

    int ploidy_ = 0;

    /**
     * @brief Index of this genotype within genotypes of the same ploidy 
     * and number of alleles.
     */
    int index_ = 0;

    /**
     * single array of alternating alleles and counts
     */
    std::vector<int> sorted_allele_counts_;

    /**
     * @brief number of different alleles in the genotype
     */
    int distinct_allele_count_ = 0;

public:
    GenotypeAlleleCounts() = default;
    GenotypeAlleleCounts(int ploidy, int index,
        std::vector<int> sorted_allele_counts, int distinct_allele_count);
    
    bool operator==(const GenotypeAlleleCounts &other) const {
        if (&other == this) return 0;
        if (ploidy_ != other.ploidy_) return false;
        if (sorted_allele_counts_.size() != other.sorted_allele_counts_.size())
            return false;
        for (size_t i = 0; i < sorted_allele_counts_.size(); ++i) {
            if (sorted_allele_counts_[i] != other.sorted_allele_counts_[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Updates the genotype allele counts to match the next genotype
     * according to the canonical ordering of PLs.
     * For example, for ploid=3 and allele_number=4, we have the following table
     * #Ploid,N_Allele  GT  GenotypeCounts
     * p3,a4   0,0,0   0,3
     * p3,a4   0,0,1   0,2,1,1
     * p3,a4   0,1,1   0,1,1,2
     * p3,a4   1,1,1   1,3
     * p3,a4   0,0,2   0,2,2,1
     * p3,a4   0,1,2   0,1,1,1,2,1
     * p3,a4   1,1,2   1,2,2,1
     * p3,a4   0,2,2   0,1,2,2
     * p3,a4   1,2,2   1,1,2,2
     * p3,a4   2,2,2   2,3
     * p3,a4   0,0,3   0,2,3,1
     * p3,a4   0,1,3   0,1,1,1,3,1
     * p3,a4   1,1,3   1,2,3,1
     * p3,a4   0,2,3   0,1,2,1,3,1
     * p3,a4   1,2,3   1,1,2,1,3,1
     * p3,a4   2,2,3   2,2,3,1
     * p3,a4   0,3,3   0,1,3,2
     * p3,a4   1,3,3   1,1,3,2
     * p3,a4   2,3,3   2,1,3,2
     * p3,a4   3,3,3   3,3
     */
    void increase();

    /**
     * @brief Updates the genotype allele counts a number of times
     * @param times number of times
     */
    void increase(int times) {
        for (int i = 0; i < times; ++i) {
            increase();
        }
    }

    int ploidy() const {
        return ploidy_;
    }

    GenotypeAlleleCounts next() const {
        GenotypeAlleleCounts this_genotype_allele_counts = *this;
        this_genotype_allele_counts.increase();
        return this_genotype_allele_counts;
    }

    int distinctAlleleCount() const {
        return distinct_allele_count_;
    }
    
    /**
     * @return double the log10 number of phased genotypes corresponding to this
     * unphased genotype.
     */
    double log10CombinationCount() {
        if (log10_combination_count_ == UNCOMPUTED_LOG_10_COMBINATION_COUNT) {
            double sum = 0.0;
            for (size_t i = 1; i < sorted_allele_counts_.size(); i+=2) {
                sum += MathUtils::log10Factorial(i);
            }
            log10_combination_count_ = MathUtils::log10Factorial(ploidy_) - sum;
        }
        return log10_combination_count_;
    }

    /**
     * @brief get the index of allele from its rank in the genotype.
     *
     * For example,
     * when there is no duplicated allele in the genotype
     * Genotype: 1/2/3, GenotypeAlleleCounts: [1, 1, 2, 1, 3, 1]
     * rank:     0/1/2
     * input rank 0 will get allele index 1,
     * input rank 1 will get allele index 2,
     * input rank 2 will get allele index 3,
     *
     * when there are duplicated alleles in the genotype
     * Genotype: 1/2/2, GenotypeAlleleCounts: [1, 1, 2, 2]
     * rank:     0/1/1
     * input rank 0 will get allele index 1,
     * input rank 1 will get allele index 2,
     * 
     * @param rank rank of the allele in the genotype
     * @return int 
     */
    int alleleIndexAt(int rank) const {
        if (rank < 0 || rank >= distinct_allele_count_) {
            std::cerr << "[GenotypeAlleleCounts::alleleIndexAt] Error! "
                << "The input rank out of range [" << 0 << ","
                << distinct_allele_count_ << ")" << std::endl;
            std::exit(1);
        }
        return sorted_allele_counts_[rank << 1];
    }

    /**
     * @brief get the rank of an allele in the genotype by its index
     * 
     * @param index allele index
     * @return int rank of index in the genotype if index is in the genotype,
     * otherwise -1 or less
     */
    int alleleRankFor(int index) const {
        if (index < 0) {
            std::cerr << "[GenotypeAlleleCounts::alleleRankFor] Error! "
                << "The input index out of range [" << 0 << ",+infinity)"
                << std::endl;
            std::exit(1);
        }
        return alleleIndexToRank(index, 0, distinct_allele_count_);
    }

    /**
     * @brief Generates a string that would represent the unphased
     * genotype with this allele counts.
     * 
     * @return std::string 
     */
    std::string toUnphasedGenotypeString() const {
        if (ploidy_ == 0) return "";
        std::vector<std::string> indices_str_vec;
        const int upper_bound = distinct_allele_count_ << 1;
        int allele_index, allele_count;
        for (int i = 0; i < upper_bound; ) {
            allele_index = sorted_allele_counts_[i++];
            allele_count = sorted_allele_counts_[i++];
            for (int j = 0; j < allele_count; ++j) {
                indices_str_vec.push_back(std::to_string(allele_index));
            }
        }
        return boost::join(indices_str_vec, VCFConstants::UNPHASED);
    }

    std::string toString() const {
        return toUnphasedGenotypeString();
    }

    /**
     * @brief get index of this genotype within all possible genotypes with 
     * this ploidy and distinct allele count
     * 
     * @return int 0 or greater
     */
    int index() const {
        return index_;
    }

    int compareTo(const GenotypeAlleleCounts &other) const {
        if (other == *this) return 0;
        if (other.ploidy_ == ploidy_) {
            return Compare(index_, other.index_);
        } else {
            return Compare(ploidy_, other.ploidy_);
        }
    }

    /**
     * @brief Returns the count of an allele in the genotype given it rank.
     * 
     * @param rank rank of the allele in the genotype
     * @return int 1 or greater
     */
    int alleleCountAt(int rank) {
        if (rank < 0 || rank >= distinct_allele_count_) {
            std::cerr << "[GenotypeAlleleCounts::alleleCountAt] Error! "
                << "The input rank out of range [" << 0 << ","
                << distinct_allele_count_ << ")" << std::endl;
            std::exit(1);
        }
        return sorted_allele_counts_[(rank << 1) + 1];
    }

    /**
     * @brief Returns the count of an allele in the genotype given it index.
     *
     * @return 0 if the allele is not present in the genotype,
     * 1 or more otherwise.
     */
    int alleleCountFor(int index) {
        const int rank = alleleRankFor(index);
        return rank < 0 ? 0 : alleleCountAt(rank);
    }

    /**
     * @brief check if this genotype contains input allele index
     * 
     * @param index input allele index
     * @return true if input allele index is contained in this genotype
     * @return false if input allele index is not contained in this genotype
     */
    bool containsAllele(int index) {
        return alleleRankFor(index) >= 0;
    }

    std::vector<int> alleleCountsByIndex(int max_allele_index) const {
        if (max_allele_index < 0)
        {
            std::cerr << "[GenotypeAlleleCounts::alleleCountsByIndex] Error!"
                "max_allele_index must >= 0" << std::endl;
            std::exit(1);
        }
        std::vector<int> result(max_allele_index+1, 0);
        copyAlleleCountsbyIndex(result, 0, 0, max_allele_index);
        return result;
    }

    /**
     * @brief Copies the sorted allele counts into an array.
     * 
     * @param dest target array
     * @param offset start offset of the target array
     */
    void copyAlleleCounts(std::vector<int> &dest, int offset) const;

    /**
     * @brief Instantiates the first genotype possible provided a total ploidy.
     * 
     * @param ploidy the ploidy of the genotype.
     * @return GenotypeAlleleCounts 
     */
    static GenotypeAlleleCounts first(const int ploidy) {
        if (ploidy < 0) {
            std::cerr << "[GenotypeAlleleCounts::first] Error!"
                "ploidy must >= 0" << std::endl;
            std::exit(1);
        }
        return ploidy == 0 ? GenotypeAlleleCounts(ploidy, 0, {}, 0) :
            GenotypeAlleleCounts(ploidy, 0, {0, ploidy}, 1);
    }

    int maximumAlleleIndex() const {
        return distinct_allele_count_ == 0 ? -1 :
            sorted_allele_counts_[(distinct_allele_count_ - 1) << 1];
    }

    int minimumAlleleIndex() const {
        return distinct_allele_count_ == 0 ? -1 : sorted_allele_counts_[0];
    }

    std::vector<int> asAlleleIndexList() const;

private:
    /**
     * @brief binary search the rank of the index in the genotype
     * 
     * @param index allele index
     * @param from first inclusive possible rank
     * @param to last exclusive possible rank
     * @return int rank of index in the genotype if index is in the genotype,
     * otherwise -1 or less. You can obtain the potential insert point
     * (within the interval [from to) ) with -result - 1
     */
    int alleleIndexToRank(int index, int from, int to) const;

    void copyAlleleCountsbyIndex(std::vector<int> &dest, int offset,
        int min_allele_index, int max_allele_index) const;
};


std::ostream &operator<<(std::ostream &out_stream,
    const GenotypeAlleleCounts &genotype_allele_counts);



#endif  // DPGT_GENOTYPE_ALLELE_COUNTS_HPP
