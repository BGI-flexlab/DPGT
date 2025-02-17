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
#ifndef VCF_READER_HPP
#define VCF_READER_HPP

#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <string>
#include <cinttypes>
#include <vector>

#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "common/simple_interval.hpp"
#include "vcf/vcf_id_table.hpp"


enum class VcfReaderQueryStatus: uint8_t {
    NOT,
    SUCCESS,
    FAIL
};


/**
 * A htslib based vcf reader for reading vcf, bgzip vcf and bcf samenessly.
 */
class VcfReader {
public:
    VcfReader() = default;
    VcfReader(const std::string &file_name, bool require_index=false) {
        Open(file_name, require_index);
    }

    VcfReader(const VcfReader &reader) = delete;
    VcfReader(VcfReader &&reader) noexcept;

    VcfReader &operator==(const VcfReader &reader) = delete;
    VcfReader &operator==(VcfReader &&reader) noexcept;

    ~VcfReader() {
        Close();
    }

    void Open(const std::string &file_name, bool require_index=false);

    void Close() {
        if (header_) bcf_hdr_destroy(header_);
        if (bcf_idx_) hts_idx_destroy(bcf_idx_);
        if (tbx_idx_) tbx_destroy(tbx_idx_);
        if (file_ptr_) hts_close(file_ptr_);
        if (itr_) hts_itr_destroy(itr_);
        if (tmps_.s) free(tmps_.s);
    }


    VcfReader &Queryi(int32_t tid, int64_t start, int64_t end);

    VcfReader &Querys(const std::string &chrom, int64_t start, int64_t end);

    VcfReader &QueryIntervals(const std::vector<SimpleInterval> &intervals);

    /**
     * Read one vcf record
     */
    bcf1_t *Read(bcf1_t *record);

    bcf_hdr_t *header() {
        return header_;
    }

    const std::vector<SimpleInterval> *getIntervals() const {
        return intervals_;
    }

    void setIdTables(VcfKeyMaps *vcf_key_maps_) {
        id_tables_ = VcfIdTables(header(), vcf_key_maps_);
    }

    VcfIdTables *getIdTables() {
        return &id_tables_;
    }

private:
    htsFile *file_ptr_ = nullptr;
    tbx_t *tbx_idx_ = nullptr;
    hts_idx_t *bcf_idx_ = nullptr;
    bcf_hdr_t *header_ = nullptr;
    hts_itr_t *itr_ = nullptr;
    bool streaming_ = true;  // streaming reader, not support random access
    kstring_t tmps_ = {0, 0, NULL};
    const std::vector<SimpleInterval> *intervals_ = nullptr; // target intervals
    std::vector<SimpleInterval>::const_iterator intervals_iter_;
    VcfIdTables id_tables_;

    VcfReaderQueryStatus query_status_ = VcfReaderQueryStatus::NOT;

    std::string getFileName() {
        if (file_ptr_ != nullptr) {
            return file_ptr_->fn;
        } else {
            return "";
        }
    }
};

#endif  // VCF_READER_HPP