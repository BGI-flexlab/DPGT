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
#include "vcf/vcf_ibuffer.hpp"
#include "htslib/vcf.h"


const int VcfIBuffer::kMaxSize = 100;

VcfIBuffer::VcfIBuffer(VcfReader *reader) {
    this->reader_ = reader;
    dual_buffer_ = (bcf1_t ***)calloc(2, sizeof(bcf1_t **));
    for (int i = 0; i < 2; ++i) {
        dual_buffer_[i] = (bcf1_t **)calloc(kMaxSize, sizeof(bcf1_t *));
        for (int j = 0; j < kMaxSize; ++j) {
            dual_buffer_[i][j] = bcf_init1();
        }
    }
    buffer_ = dual_buffer_[0];
    FillBuffer();
    buffer_idx_ = 0;
}

VcfIBuffer::VcfIBuffer(VcfIBuffer &&vcf_ibuffer) noexcept {
    this->reader_ = vcf_ibuffer.reader_;
    this->cur_ = vcf_ibuffer.cur_;
    this->buffer_ = vcf_ibuffer.buffer_;
    this->dual_buffer_ = vcf_ibuffer.dual_buffer_;
    this->buffer_idx_ = vcf_ibuffer.buffer_idx_;

    vcf_ibuffer.reader_ = nullptr;
    vcf_ibuffer.buffer_ = nullptr;
    vcf_ibuffer.dual_buffer_ = nullptr;
}

VcfIBuffer & VcfIBuffer::operator==(VcfIBuffer &&vcf_ibuffer) noexcept {
    if (&vcf_ibuffer != this) {
        Free();
        
        this->reader_ = vcf_ibuffer.reader_;
        this->cur_ = vcf_ibuffer.cur_;
        this->buffer_ = vcf_ibuffer.buffer_;
        this->dual_buffer_ = vcf_ibuffer.dual_buffer_;
        this->buffer_idx_ = vcf_ibuffer.buffer_idx_;

        vcf_ibuffer.reader_ = nullptr;
        vcf_ibuffer.buffer_ = nullptr;
        vcf_ibuffer.dual_buffer_ = nullptr;
    }
    return *this;
}

VcfIBuffer::~VcfIBuffer() {
    Free();
}

void VcfIBuffer::FillBuffer() {
    cur_ = 0;
    for (int i = 0; i < kMaxSize; ++i) {
        if (!(reader_->Read(buffer_[i]))) {
            bcf_destroy1(buffer_[i]);
            buffer_[i] = nullptr;
            break;
        }
    }
}

void VcfIBuffer::ReinitializeBuffer() {
    cur_ = 0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < kMaxSize; ++j) {
            if (dual_buffer_[i][j] != NULL) {
                bcf_clear(dual_buffer_[i][j]);
            } else {
                dual_buffer_[i][j] = bcf_init1();
            }
        }
    }
    buffer_ = dual_buffer_[0];
    FillBuffer();
    buffer_idx_ = 0;
}

void VcfIBuffer::Queryi(int32_t tid, int64_t start, int64_t end) {
    reader_->Queryi(tid, start, end);
    ReinitializeBuffer();
}

void VcfIBuffer::Querys(const std::string &chrom, int64_t start, int64_t end) {
    reader_->Querys(chrom, start, end);
    ReinitializeBuffer();
}

void VcfIBuffer::QueryIntervals(
    const std::vector<SimpleInterval> &intervals)
{
    reader_->QueryIntervals(intervals);
    ReinitializeBuffer();
}

bool VcfIBufferLess::operator()(VcfIBuffer *left, VcfIBuffer *right) const {
    auto func = [](bcf1_t *a, bcf1_t *b) -> int {
        int cmp = 0;
        if ((cmp = Compare(a->rid, b->rid)) != 0) {
            return cmp;
        } else if ((cmp = Compare(a->pos, b->pos)) != 0) {
            return cmp;
        } else {
            return 0;
        }
    };

    if (left->First() != nullptr && right->First() != nullptr) {
        return func(left->First(), right->First()) < 0;
    } else if (left->First() == nullptr && right->First() != nullptr) {
        return false;
    } else if (left->First() != nullptr && right->First() == nullptr) {
        return true;
    } else {
        return false;
    }
}


bool VcfIBufferGreater::operator()(VcfIBuffer *left, VcfIBuffer *right) const {
    auto func = [](bcf1_t *a, bcf1_t *b) -> int {
        int cmp = 0;
        if ((cmp = Compare(a->rid, b->rid)) != 0) {
            return cmp;
        } else if ((cmp = Compare(a->pos, b->pos)) != 0) {
            return cmp;
        } else {
            return 0;
        }
    };

    if (left->First() != nullptr && right->First() != nullptr) {
        return func(left->First(), right->First()) > 0;
    } else if (left->First() == nullptr && right->First() != nullptr) {
        return true;
    } else if (left->First() != nullptr && right->First() == nullptr) {
        return false;
    } else {
        return false;
    }
}
