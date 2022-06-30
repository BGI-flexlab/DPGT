#ifndef DPGT_VCF_ID_TABLE_HPP
#define DPGT_VCF_ID_TABLE_HPP


#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>
#include "htslib/vcf.h"


class VcfKeyMap {
public:
    std::vector<std::string> keys;
    std::unordered_map<std::string, int> key_indices;

    virtual ~VcfKeyMap() {}

    size_t size() const {
        return keys.size();
    }

    std::string getKeyByIdx(int i) const
    {
        return keys[i];
    }

    int getIdxByKey(const char *k) const {
        auto itr = key_indices.find(k);
        if (itr != key_indices.end()) {
            return itr->second;
        }
        return -1;
    }
};


class VcfSampleKeyMap: public VcfKeyMap {
public:
    VcfSampleKeyMap() = default;
    VcfSampleKeyMap(bcf_hdr_t *hdr);

    virtual ~VcfSampleKeyMap() override {}
};


class VcfFormatKeyMap: public VcfKeyMap {
public:
    int AD_IDX = -1;
    int DP_IDX = -1;
    int GQ_IDX = -1;
    int PL_IDX = -1;
    int MIN_DP_IDX = -1;

    VcfFormatKeyMap() = default;
    VcfFormatKeyMap(bcf_hdr_t *hdr);

    virtual ~VcfFormatKeyMap() override {}
};


class VcfIdTable {
public:
    std::vector<int> table;
    VcfKeyMap *key_map;

    virtual ~VcfIdTable() {}
};


class VcfSampleIdTable: public VcfIdTable {
public:
    VcfSampleIdTable() = default;
    VcfSampleIdTable(bcf_hdr_t *hdr, VcfSampleKeyMap *key_map_);

    virtual ~VcfSampleIdTable() override {}
};


class VcfFormatIdTable: public VcfIdTable {
public:
    VcfFormatIdTable() = default;
    VcfFormatIdTable(bcf_hdr_t *hdr, VcfFormatKeyMap *key_map_);

    virtual ~VcfFormatIdTable() override {}
};


class VcfKeyMaps {
public:
    VcfSampleKeyMap sample_key_map;
    VcfFormatKeyMap format_key_map;
    VcfKeyMaps() = default;
    VcfKeyMaps(bcf_hdr_t *hdr) {
        sample_key_map = VcfSampleKeyMap(hdr);
        format_key_map = VcfFormatKeyMap(hdr);
    }
};


class VcfIdTables {
public:
    VcfSampleIdTable sample_id_table;
    VcfFormatIdTable format_id_table;

    VcfIdTables() = default;

    VcfIdTables(bcf_hdr_t *hdr, VcfKeyMaps *key_maps)
    {
        sample_id_table = VcfSampleIdTable(hdr, &key_maps->sample_key_map);
        format_id_table = VcfFormatIdTable(hdr, &key_maps->format_key_map);
    }
};


#endif // DPGT_VCF_ID_TABLE_HPP
