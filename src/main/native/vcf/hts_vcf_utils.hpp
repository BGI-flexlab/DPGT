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
#ifndef DPGT_HTS_VCF_UTILS_HPP
#define DPGT_HTS_VCF_UTILS_HPP

#include <string>
#include <utility>
#include <vector>
#include "htslib/vcf.h"


/**
 * @brief rename bcf header samples
 * 
 * @param header old header
 * @param old_idx_new_name_pairs old index and new sample name pairs
 * @return bcf_hdr_t* new header
 */
bcf_hdr_t *bcf_hdr_rename_samples(bcf_hdr_t *header,
    const std::vector<std::pair<int, std::string>> &old_idx_new_name_pairs);

/**
 * @brief rename duplicated samples by add suffix(-int) to them
 * 
 * @param headers bcf headers to opperate on
 */
void bcf_hdr_rename_dup_samples(std::vector<bcf_hdr_t *> headers,
    bool allow_dup_samples = true);


/**
 * @brief merge vcf headers and add samples
 *
 * assume that input headers do not have duplicated samples(
 * duplicated samples can be renamed by bcf_hdr_rename_dup_samples)
 * 
 * @return bcf_hdr_t* merged vcf header, need to be freed after use.
 */
bcf_hdr_t *bcf_hdr_merge_add_samples(const std::vector<bcf_hdr_t *> headers);


#endif  // DPGT_HTS_VCF_UTILS_HPP
