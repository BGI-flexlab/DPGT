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
#include "genotype_allele_counts.hpp"
#include <cstddef>
#include <cstdlib>
#include <vector>
#include <iostream>


 const double GenotypeAlleleCounts::UNCOMPUTED_LOG_10_COMBINATION_COUNT = -1;

template <typename T>
static void copyVector(const std::vector<T> &src, int src_pos,
    std::vector<T> &dest, int dest_pos, int length)
{
    if (src_pos < 0) {
        std::cerr << "[copyVector] out of range error! src_pos < 0"
            << std::endl;
        std::exit(1);
    }
    if (dest_pos < 0) {
        std::cerr << "[copyVector] out of range error! dest_pos < 0"
            << std::endl;
        std::exit(1);
    }
    const int to_copy_src_size = length + src_pos;
    if (to_copy_src_size > (int)src.size())
    {
        std::cerr << "[copyVector] out of range error! "
            "length + src_pos > src size" << std::endl;
        std::exit(1);
    }
    const int new_dest_size = length + dest_pos;
    if (new_dest_size > (int)dest.size())
    {
        std::cerr << "[copyVector] out of range error! "
            "length + src_pos > src size" << std::endl;
        std::exit(1);
    }

    for (int i = 0; i < length; ++i) {
        dest[dest_pos++] = src[src_pos++];
    }
}


GenotypeAlleleCounts::GenotypeAlleleCounts(int ploidy, int index,
    std::vector<int> sorted_allele_counts, int distinct_allele_count):
    ploidy_(ploidy), index_(index),
    sorted_allele_counts_(std::move(sorted_allele_counts)),
    distinct_allele_count_(distinct_allele_count) {}


void GenotypeAlleleCounts::increase() {
    // if the ploidy is zero there is only one possible genotype.
    if (distinct_allele_count_ == 0) return;

    if (distinct_allele_count_ == 1) {
        if (ploidy_ == 1) {
            ++sorted_allele_counts_[0];
        } else {
            if (sorted_allele_counts_.size() < 4) {
                sorted_allele_counts_.resize(4);
            }
            sorted_allele_counts_[2] = sorted_allele_counts_[0] + 1;
            sorted_allele_counts_[3] = 1;
            sorted_allele_counts_[0] = 0;
            sorted_allele_counts_[1] = ploidy_ - 1;
            distinct_allele_count_ = 2;
        }
    } else {
        const int allele0 = sorted_allele_counts_[0];
        const int freq0 = sorted_allele_counts_[1];
        const int allele1 = sorted_allele_counts_[2];
        const int allele0_plus1 = allele0 + 1;
        const bool allele0_1_are_consecutive = allele0_plus1 == allele1;
        const std::vector<int> new_sorted_allele_counts_;
        // The rest of the sorted allele counts array contains junk
        const int sorted_allele_counts_length = distinct_allele_count_ << 1;


        if (freq0 == 1) {   // in this case allele0 wont be present in the result and all is frequency should go to allele0 + 1.
            if (allele0_1_are_consecutive) {  // need just to remove the first allele and add 1 to the frequency of the second (freq1 += 1).
                sorted_allele_counts_.erase(sorted_allele_counts_.begin(), sorted_allele_counts_.begin()+2); // shift left the first component away.
                sorted_allele_counts_[1]++; // freq1 has become freq0.
                distinct_allele_count_--;
            } else  // just need to mutate allele0 to allele0 + 1.
            {
                sorted_allele_counts_[0] = allele0_plus1;
            }
        } else { // && freq0 > 1 as per sorted_allele_counts_ format restrictions. In this case allele0 will mutated to '0' with frequency decreased by 1.
            if (allele0_1_are_consecutive) { // we don't need to add a component for allele0 + 1 since it already exists.
                sorted_allele_counts_[0] = 0;
                sorted_allele_counts_[1] = freq0 - 1;
                sorted_allele_counts_[3]++;
            } else { // we need to insert allele0 + 1 in the sorted-allele-counts array and give it frequency 1.
                if (sorted_allele_counts_.size() < (size_t)sorted_allele_counts_length + 2) // make room for the new component.
                {
                    sorted_allele_counts_.resize(sorted_allele_counts_length + 2);
                }
                copyVector(sorted_allele_counts_, 2, sorted_allele_counts_, 4,
                    sorted_allele_counts_length - 2);
                
                sorted_allele_counts_[0] = 0;
                sorted_allele_counts_[1] = freq0 - 1;
                sorted_allele_counts_[2] = allele0_plus1;
                sorted_allele_counts_[3] = 1;
                distinct_allele_count_++;
            }
        }
    }

    ++index_;
    log10_combination_count_ = UNCOMPUTED_LOG_10_COMBINATION_COUNT;
}


int GenotypeAlleleCounts::alleleIndexToRank(int index, int from, int to) const {
    int mid;
    int mid_index;
    while (from < to) {
        mid = (to + from) >> 1;
        mid_index = sorted_allele_counts_[mid << 1];
        if (mid_index == index) return mid;
        if (mid_index < index) {
            from = mid + 1;
        } else {
            to = mid;
        }
    }

    return -to - 1;

// not use recursive implementation
#if 0    
    if (from == to - 1) {
        const int only_index = sorted_allele_counts_[from << 1];
        return only_index == index ? from :
            (only_index > index ? -from - 1 : -to - 1);
    }

    const int mid = (to + from) >> 1;
    const int mid_index = sorted_allele_counts_[mid << 1];
    if (mid_index == index) {
        return mid;
    } else if (mid_index < index) {
        return alleleIndexToRank(index, mid + 1, to);
    } else {
        return alleleIndexToRank(index, 0, mid);
    }
#endif
}


void GenotypeAlleleCounts::copyAlleleCountsbyIndex(
    std::vector<int> &dest, int offset,
    int min_allele_index, int max_allele_index) const
{
    const int min_allele_rank = alleleRankFor(min_allele_index);
    const int max_allele_rank = alleleRankFor(max_allele_index);

    const int start_rank =
        min_allele_rank < 0 ? (-min_allele_rank - 1) : min_allele_rank;
    const int end_rank =
        max_allele_rank < 0 ? (-max_allele_rank - 2) : max_allele_rank;

    int next_index = min_allele_index;
    int next_rank = start_rank;
    int next_sorted_allele_counts_offset = next_rank << 1;
    int next_dest_offset = offset;
    while (next_rank++ <= end_rank) {
        const int allele_index =
            sorted_allele_counts_[next_sorted_allele_counts_offset++];
        while (allele_index > next_index) {
            dest[next_dest_offset++] = 0;
            next_index++;
        }
        dest[next_dest_offset++] =
            sorted_allele_counts_[next_sorted_allele_counts_offset++];
        next_index++;
    }
    while (next_index++ <= max_allele_index) {
        dest[next_dest_offset++] = 0;
    }
}

void GenotypeAlleleCounts::copyAlleleCounts(
    std::vector<int> &dest, int offset) const
{
    if (offset < 0) {
        std::cerr << "[GenotypeAlleleCounts::copyAlleleCounts] Error!"
            "offset must >= 0" << std::endl;
        std::exit(1);
    }

    if (offset + sorted_allele_counts_.size() > dest.size()) {
        std::cerr << "[GenotypeAlleleCounts::copyAlleleCounts] Error!"
            "input vector do not have enough capacity" << std::endl;
        std::exit(1);
    }

    copyVector(sorted_allele_counts_, 0, dest, offset,
        sorted_allele_counts_.size());
}

std::vector<int> GenotypeAlleleCounts::asAlleleIndexList() const {
    std::vector<int> result(ploidy_, 0);
    int offset = 0;
    for (int i = 0; i < distinct_allele_count_; ++i) {
        const int allele_index = sorted_allele_counts_[i];
        const int allele_count = sorted_allele_counts_[i+i];
        for (int j = 0; j < allele_count; ++j) {
            result[offset++] = allele_index;
        }
    }
    return result;
}

std::ostream &operator<<(std::ostream &out_stream,
    const GenotypeAlleleCounts &genotype_allele_counts)
{
    out_stream << genotype_allele_counts.toString();
    return out_stream;
}
