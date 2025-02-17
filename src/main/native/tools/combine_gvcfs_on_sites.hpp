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
#ifndef DPGT_COMBINE_GVCFS_ON_SITES_HPP
#define DPGT_COMBINE_GVCFS_ON_SITES_HPP

#include <cstdint>
#include <list>
#include <memory>
#include <string>
#include <vector>
#include "common/cached_ref_seq.hpp"
#include "vcf/multi_vcf_reader.hpp"
#include "vcf/variant_context.hpp"
#include "vcf/variant_context_utils.hpp"
#include "variant/variant_site_set.hpp"
#include "vcf/variant_builder.hpp"
#include "vcf/vcf_writer.hpp"
#include "tools/reference_confident_variant_merger.hpp"
#include "vcf/gatk_vcf_constants.hpp"

/**
 * @brief Combine gvcfs on input sites
 */
class CombineGVCFsOnSites {
public:
    CombineGVCFsOnSites(MultiVcfReader *reader,
        CachedRefSeq *cached_ref_seq, VariantSiteSet *vs_set,
        std::string chrom, VcfWriter *vcf_writer):
        reader_(reader), merged_header_(reader->header()),
        cached_ref_seq_(cached_ref_seq), vs_set_(vs_set),
        chrom_(std::move(chrom)), vcf_writer_(vcf_writer)
    {
        bcf_hdr_append(merged_header_,
            GATKVCFConstants::RAW_MAPPING_QUALITY_WITH_DEPTH_KEY_LINE.c_str());
        cur_site_ = vs_set->firstSite();
        tid_ = bcf_hdr_name2id(getMergedHeader(), chrom_.c_str());
        merged_samples_ = createMergedSamples();
        vcf_writer_->writeHeader(getMergedHeader());
    }

    void apply(std::shared_ptr<VariantContext> &&var);
    void apply(int64_t site,
        std::list<std::shared_ptr<VariantContext>> &variants);
    void finalize();
    void run();

private:
    MultiVcfReader *reader_ = nullptr;
    bcf_hdr_t *merged_header_ = nullptr;
    CachedRefSeq *cached_ref_seq_ = nullptr;
    VariantSiteSet *vs_set_ = nullptr;
    std::string chrom_;
    int32_t tid_;
    int64_t cur_site_ = -1;
    std::list<std::shared_ptr<VariantContext>> variants_overlap_cur_site_;
    kstring_t tmp_var_ks_ = {0, 0, NULL};
    VcfWriter *vcf_writer_ = nullptr;

    ReferenceConfidentVariantMerger reference_confident_var_merger_;

    std::vector<std::string> merged_samples_;
    // merged header related functions
    bcf_hdr_t *getMergedHeader() {
        return merged_header_;
    }

    std::vector<std::string> createMergedSamples() {
        bcf_hdr_t *header = getMergedHeader();
        std::vector<std::string> result;
        for (int i = 0; i < bcf_hdr_nsamples(header); ++i) {
            result.push_back(header->samples[i]);
        }
        return result;
    }

    VcfKeyMaps *getKeyMaps() {
        return reader_->KeyMaps();
    }

    std::vector<std::string> getMergedSamples() {
        return merged_samples_;
    }

    int getMergedSamplesNumber() {
        return merged_samples_.size();
    }

    bool containsTrueAltAllele(
        const std::list<std::shared_ptr<VariantContext>> &variants)
    {
        for (auto const &vc: variants) {
            if (vc->getNAllele() > 2) return true;
        }
        return false;
    }

    VariantStringWithPos referenceBlockMerge(
        std::vector<std::shared_ptr<VariantContext>> &variants, int64_t site);
};



#endif  // DPGT_COMBINE_GVCFS_ON_SITES_HPP
