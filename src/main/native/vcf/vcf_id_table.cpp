#include <algorithm>
#include "vcf/vcf_id_table.hpp"


VcfSampleKeyMap::VcfSampleKeyMap(bcf_hdr_t *hdr): VcfKeyMap()
{
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        keys.push_back(hdr->samples[i]);
    }
    key_indices.reserve(keys.size());
    for (int i = 0; i < static_cast<int>(keys.size()); ++i) {
        key_indices[keys[i]] = i;
    }
}


VcfFormatKeyMap::VcfFormatKeyMap(bcf_hdr_t *hdr): VcfKeyMap() {
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


VcfSampleIdTable::VcfSampleIdTable(bcf_hdr_t *hdr,
    VcfSampleKeyMap *key_map_): VcfIdTable()
{
    key_map = key_map_;
    std::vector<int> indices1;
    std::vector<int> indices2;
    int max = 0;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        auto itr = key_map->key_indices.find(hdr->samples[i]);
        if (itr != key_map->key_indices.end()) {
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


VcfFormatIdTable::VcfFormatIdTable(bcf_hdr_t *hdr,
    VcfFormatKeyMap *key_map_): VcfIdTable()
{
    key_map = key_map_;
    std::vector<int> indices1;
    std::vector<int> indices2;
    int max = 0;
    for (int32_t i = 0; i < hdr->n[BCF_DT_ID]; ++i) {
        const char *k = hdr->id[BCF_DT_ID][i].key;
        auto itr = key_map->key_indices.find(k);
        if (itr != key_map->key_indices.end()) {
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
