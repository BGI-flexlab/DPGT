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
#include <algorithm>
#include <cstddef>
#include <memory>
#include <sys/types.h>
#include <vector>
#include "vcf/hts_vcf_utils.hpp"
#include "vcf/multi_vcf_reader.hpp"
#include "htslib/vcf.h"
#include "vcf/variant_context.hpp"


MultiVcfReader::MultiVcfReader(std::vector<std::string> files,
    bool require_index):files_(std::move(files))
{
    for (size_t i = 0; i < files_.size(); ++i) {
        VcfReader *reader = new VcfReader(files_[i], require_index);
        readers_.push_back(reader);
        VcfIBuffer *buffer = new VcfIBuffer(reader);
        if (buffer->First() != nullptr) {
            buffers_.push(buffer);
            buffers_vec_.push_back(buffer);
        } else {
            delete buffer;
        }
    }

    std::vector<bcf_hdr_t *> headers;
    headers.reserve(readers_.size());
    for (auto &r: readers_) headers.push_back(r->header());
    header_ = bcf_hdr_merge_add_samples(headers);

    // format key map for merged header
    vcf_key_maps_ = new VcfKeyMaps(header_);

    for (auto &b: buffers_vec_) {
        b->setIdTables(vcf_key_maps_);
    }

}


MultiVcfReader::MultiVcfReader(MultiVcfReader &&multi_vcf_reader) noexcept {
    this->files_ = std::move(multi_vcf_reader.files_);
    this->readers_ = std::move(multi_vcf_reader.readers_);
    this->buffers_ = std::move(multi_vcf_reader.buffers_);
    this->buffers_vec_ = std::move(multi_vcf_reader.buffers_vec_);
}

MultiVcfReader &MultiVcfReader::operator==(
    MultiVcfReader &&multi_vcf_reader) noexcept
{
    if (&multi_vcf_reader != this) {
        Free();
        this->files_ = std::move(multi_vcf_reader.files_);
        this->readers_ = std::move(multi_vcf_reader.readers_);
        this->buffers_ = std::move(multi_vcf_reader.buffers_);
        this->buffers_vec_ = std::move(multi_vcf_reader.buffers_vec_);
    }
    return *this;
}

MultiVcfReader::~MultiVcfReader() {
    Free();
}

void MultiVcfReader::ReloadBuffers() {
    std::vector<VcfIBuffer *> tmp;
    tmp.reserve(buffers_.size());
    while (!buffers_.empty()) {
        VcfIBuffer *b = buffers_.top();
        if (b->First() != nullptr) tmp.push_back(b);
        buffers_.pop();
    }
    for (VcfIBuffer *b: tmp) buffers_.push(b);
}

std::shared_ptr<VariantContext> MultiVcfReader::First() {
    if (buffers_.empty()) {
        return nullptr;
    }
    bcf1_t *record = buffers_.top()->First();
    bcf_hdr_t *header = buffers_.top()->header();
    VcfIdTables *id_table = buffers_.top()->getIdTables();
    auto vc = new VariantContext(record, header, id_table);
    return std::shared_ptr<VariantContext>(vc);
}

std::shared_ptr<VariantContext> MultiVcfReader::Read() {
    if (buffers_.empty()) {
        return nullptr;
    }
    VcfIBuffer *b = buffers_.top();
    bcf_hdr_t *header = buffers_.top()->header();
    VcfIdTables *id_table = buffers_.top()->getIdTables();
    bcf1_t *record = b->Read();

    buffers_.pop();
    if (b->First() != nullptr) {
        buffers_.push(b);
    }

    auto vc = new VariantContext(record, header, id_table);
    return std::shared_ptr<VariantContext>(vc);
}

bool MultiVcfReader::IsEOF() {
    return First() == nullptr;
}


void MultiVcfReader::Queryi(int32_t tid, int64_t start, int64_t end) {
    for (VcfIBuffer *r: buffers_vec_) r->Queryi(tid, start, end);
    ReloadBuffers();
}

void MultiVcfReader::Querys(const std::string &chrom,
    int64_t start, int64_t end)
{
    for (VcfIBuffer *r: buffers_vec_) r->Querys(chrom, start, end);
    ReloadBuffers();
}

void MultiVcfReader::QueryIntervals(
    const std::vector<SimpleInterval> &intervals)
{
    for (VcfIBuffer *r: buffers_vec_) r->QueryIntervals(intervals);
    ReloadBuffers();
}

bcf_hdr_t *MultiVcfReader::header() {
    return header_;
}
