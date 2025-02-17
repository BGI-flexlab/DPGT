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
#ifndef DPGT_VCF_HEADER_ID_TABLE_HPP
#define DPGT_VCF_HEADER_ID_TABLE_HPP

#include "htslib/vcf.h"
#include "htslib/khash.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>


class VcfFmtKeyMap {
public:
    std::vector<std::string> keys;
    std::unordered_map<std::string, int> key_indices;
    int AD_IDX = -1;
    int DP_IDX = -1;
    int GQ_IDX = -1;
    int PL_IDX = -1;
    int MIN_DP_IDX = -1;

    VcfFmtKeyMap() = default;
    VcfFmtKeyMap(bcf_hdr_t *hdr);

    size_t size() const {
        return keys.size();
    }

    std::string getKeyByIdx(int i) const;

    int getIdxByKey(const char *k) const;
};

class VcfHeaderFmtIdTable {
public:
    // mapping from individual vcf header FORMAT key id to merged vcf header 
    // FORMAT key index(not the same as merged vcf header key id)
    std::vector<int> table;

    VcfFmtKeyMap *merged_vcf_fmt_key_map;

    VcfHeaderFmtIdTable() = default;
    VcfHeaderFmtIdTable(bcf_hdr_t *hdr, VcfFmtKeyMap *key_map);
};


#endif  // DPGT_VCF_HEADER_ID_TABLE_HPP
