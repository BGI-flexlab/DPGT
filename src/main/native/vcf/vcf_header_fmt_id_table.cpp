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
#include "vcf/vcf_header_fmt_id_table.hpp"
#include <string>


std::string VcfFmtKeyMap::getKeyByIdx(int i) const {
    return keys[i];
}

int VcfFmtKeyMap::getIdxByKey(const char *k) const
{
    auto itr = key_indices.find(k);
    if (itr != key_indices.end()) {
        return itr->second;
    }
    return -1;
}

VcfFmtKeyMap::VcfFmtKeyMap(bcf_hdr_t *hdr) {
    for (int i = 0; i < hdr->nhrec; ++i) {
        if (hdr->hrec[i]->type != BCF_HL_FMT) continue;
        for (int j = 0; j < hdr->hrec[i]->nkeys; ++j) {
            if (strcmp(hdr->hrec[i]->keys[j], "ID") == 0) {
                keys.push_back(hdr->hrec[i]->vals[j]);
                break;
            }
        }
    }
    std::sort(keys.begin(), keys.end());
    key_indices.reserve(keys.size());
    for (int i = 0; i < static_cast<int>(keys.size()); ++i) {
        key_indices[keys[i]] = i;
    }

    AD_IDX = getIdxByKey("AD");
    DP_IDX = getIdxByKey("DP");
    GQ_IDX = getIdxByKey("GQ");
    PL_IDX = getIdxByKey("PL");
    MIN_DP_IDX = getIdxByKey("MIN_DP");
}

VcfHeaderFmtIdTable::VcfHeaderFmtIdTable(bcf_hdr_t *hdr,
    VcfFmtKeyMap *key_map)
{
    merged_vcf_fmt_key_map = key_map;
    std::vector<int> indices1;
    std::vector<int> indices2;
    int max = 0;
    for (int32_t i = 0; i < hdr->n[BCF_DT_ID]; ++i) {
        const char *k = hdr->id[BCF_DT_ID][i].key;
        auto itr = merged_vcf_fmt_key_map->key_indices.find(k);
        if (itr != merged_vcf_fmt_key_map->key_indices.end()) {
            indices1.push_back(i);
            indices2.push_back(itr->second);
            if (i > max) max = i;
        }
    }

    table = std::vector<int>(max+1, -1);
    for (size_t i = 0; i < indices1.size(); ++i) {
        table[indices1[i]] = indices2[i];
    }
}
