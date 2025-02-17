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
#ifndef VCF_IBUFFER_HPP
#define VCF_IBUFFER_HPP

#include <cstdlib>
#include "htslib/hts_endian.h"
#include "htslib/vcf.h"
#include "vcf/vcf_reader.hpp"
#include "vcf/vcf_id_table.hpp"
#include "common/math_utils.hpp"


/**
 * vcf input buffer.
 */
class VcfIBuffer {
public:
    static const int kMaxSize;

    VcfIBuffer(VcfReader *reader);

    VcfIBuffer(const VcfIBuffer &vcf_ibuffer) = delete;
    VcfIBuffer(VcfIBuffer &&vcf_ibuffer) noexcept;

    VcfIBuffer &operator==(const VcfIBuffer &vcf_ibuffer) = delete;
    VcfIBuffer &operator==(VcfIBuffer &&vcf_ibuffer) noexcept;

    ~VcfIBuffer();

    bcf1_t *First() {
        if (cur_ >= kMaxSize) {
            buffer_idx_ = 1 - buffer_idx_;
            buffer_ = dual_buffer_[buffer_idx_];
            FillBuffer();
        }
        return buffer_[cur_];
    }

    bcf1_t *Read() {
        if (cur_ >= kMaxSize) {
            buffer_idx_ = 1 - buffer_idx_;
            buffer_ = dual_buffer_[buffer_idx_];
            FillBuffer();
        }
        return buffer_[cur_++];
    }

    bool IsEOF() {
        return First() == nullptr;
    }

    void Queryi(int32_t tid, int64_t start, int64_t end);

    void Querys(const std::string &chrom, int64_t start, int64_t end);

    void QueryIntervals(const std::vector<SimpleInterval> &intervals);

    VcfReader *reader() const {
        return reader_;
    }

    bcf_hdr_t *header() const {
        if (reader_ == nullptr) return nullptr;
        return reader_->header();
    }

    const std::vector<SimpleInterval> *getIntervals() const {
        if (reader_ == nullptr) return nullptr;
        return reader_->getIntervals();
    }

    void setIdTables(VcfKeyMaps *vcf_key_maps_) {
        id_tables_ = VcfIdTables(header(), vcf_key_maps_);
    }

    VcfIdTables *getIdTables() {
        return &id_tables_;
    }

private:
    VcfReader *reader_ = nullptr;
    int cur_ = 0;
    bcf1_t **buffer_ = nullptr;         // pointer to dual_buffer_ element
    bcf1_t ***dual_buffer_ = nullptr;   // array of 2 bcf1_t **
    uint32_u buffer_idx_ = 0;           // a switch
    VcfIdTables id_tables_;

    void FillBuffer();

    void ReinitializeBuffer();

    void Free() {
        if (dual_buffer_) {
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < kMaxSize; ++j) {
                bcf_destroy1(dual_buffer_[i][j]);
            }
            free(dual_buffer_[i]);
        }
            free(dual_buffer_);
        }
    }
};


struct VcfIBufferLess {
    bool operator()(VcfIBuffer *left, VcfIBuffer *right) const;
};


struct VcfIBufferGreater {
    bool operator()(VcfIBuffer *left, VcfIBuffer *right) const;
};


#endif  // VCF_IBUFFER_HPP
