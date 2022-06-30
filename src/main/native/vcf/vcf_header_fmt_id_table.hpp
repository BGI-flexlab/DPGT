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
