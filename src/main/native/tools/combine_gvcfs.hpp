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
#ifndef DPGT_COMBINE_GVCFS_HPP
#define DPGT_COMBINE_GVCFS_HPP

#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <list>
#include <set>
#include "common/simple_interval.hpp"
#include "common/utils.hpp"
#include "htslib/kstring.h"
#include "tools/reference_confident_variant_merger.hpp"
#include "vcf/variant_context.hpp"
#include "vcf/variant_context_utils.hpp"
#include "common/reference_context.hpp"
#include "tools/multi_variant_group_on_start.hpp"
#include "boost/dynamic_bitset.hpp"
#include "vcf/vcf_writer.hpp"


class RelativeBitSet {
private:
    boost::dynamic_bitset<> data_;
    int64_t start_ = 0;
    int64_t size_;
public:
    RelativeBitSet() = default;
    RelativeBitSet(int64_t start, int64_t end): start_(start) {
        size_ = end - start + 1;
        data_ = boost::dynamic_bitset<>(size_);
    }

    RelativeBitSet &set(int64_t val) {
        Utils::validateArg(val >= start_,
            "input value can not less than start");
        data_.set(val - start_);
        return *this;
    }

    RelativeBitSet &reset(int64_t val) {
        Utils::validateArg(val >= start_,
            "input value can not less than start");
        data_.reset(val - start_);
        return *this;
    }

    bool test(int64_t val) {
        Utils::validateArg(val >= start_,
            "input value can not less than start");
        return data_.test(val - start_);
    }

    bool testIndex(int64_t i) {
        return data_.test(i);
    }

    int64_t get(int64_t i) {
        return i + start_;
    }

    int64_t size() const {
        return size_;
    }
};


/**
 * @brief a tool for merge gvcfs
 */
class CombineGVCFs: public MultiVariantGroupOnStart {
private:
    static ReferenceConfidentVariantMerger reference_confident_var_merger_;
    std::list<std::shared_ptr<VariantContext>> variants_overlap_current_merge_;
    boost::dynamic_bitset<> sample_indices_;
    SimpleInterval prev_pos_;
    char ref_after_prev_pos_;
    ReferenceContext stored_reference_context_;
    VcfWriter *vcf_writer_ = nullptr;
    kstring_t tmp_var_ks_ = {0, 0, NULL};

    /**
     * @brief If > 0, reference bands will be broken up at genomic positions
     * that are multiples of this number
     * For example, a value of 10,000 would mean that we would ensure that
     * no bands(gvcf block) span across chr1:10000, chr1:20000, etc.
     * default value 0 means that not break bands
     */
    int64_t multiple_at_witch_to_break_bands_ = 0;

    /**
     * @brief if true, output gvcf block at every position
     */
    bool use_bp_resolution_ = false;

    template<typename T>
    static int64_t maxVariantContextsEnd(const T &begin, const T &end) {
        int64_t max = 0;
        for (T itr = begin; itr != end; ++itr) {
            if ((*itr)->getEnd() > max) max = (*itr)->getEnd();
        }
        return max;
    }

    void createIntermediateVariants(const SimpleInterval &interval_to_close);

    /**
     * @brief Get the Intermediate Stop Sites by break at multiples of 
     * break_bands_multiple
     * For example, given interval 1:25-48 and break_bands_multiple 4,
     * output breaks: 27, 31, 35, 39, 43, 47
     * the start of interval to close is not included in output breaks, but end
     * of interval of closed is included.
     * 
     * @param interval_to_close 
     * @param break_bands_multiple 
     * @return std::set<int64_t> 
     */
    RelativeBitSet getIntermediateStopSites(
        const SimpleInterval &interval_to_close,
        int64_t break_bands_multiple);
    
    /**
     * @brief resize stored_reference_context to cover at least as much as
     * interval_to_close
     */
    void resizeReferenceIfNeeded(const SimpleInterval &interval_to_close) {
        int64_t left = stored_reference_context_.getStart() -
            interval_to_close.start;
        int64_t right = interval_to_close.end -
            stored_reference_context_.getEnd();
        stored_reference_context_.setWindow(
            std::max(1L, left), std::max(1L, right));
    }

    void mergeWithNewVCs(
        const std::vector<std::shared_ptr<VariantContext>> &variants,
        ReferenceContext &reference_context);

    bool okayToSkipThisSite(
        const std::vector<std::shared_ptr<VariantContext>> &variants,
        const ReferenceContext &reference_context);
    
    static std::set<std::string> getSamples(
        const std::vector<std::shared_ptr<VariantContext>> &variants);
    
    static boost::dynamic_bitset<> getSampleIndices(
        const std::vector<std::shared_ptr<VariantContext>> &variants);

    void endPreviousStates(const SimpleInterval &pos,
        const std::string &ref_bases,
        const std::vector<std::shared_ptr<VariantContext>> &variants,
        bool force_output_at_current_pos);
    
    /**
     * @param variants input variants
     * @return true if there are true(exclude <NON-REF>) alternative alleles in
     * input variants
     */
    bool containsTrueAltAllele(
        const std::vector<std::shared_ptr<VariantContext>> &variants);

    VariantStringWithPos referenceBlockMerge(
        std::vector<std::shared_ptr<VariantContext>> &variants, int64_t end);

public:
    CombineGVCFs(MultiVcfReader *reader,
        CachedRefSeq *cached_ref_seq, VcfWriter *vcf_writer):
        MultiVariantGroupOnStart(reader, cached_ref_seq),
        vcf_writer_(vcf_writer)
    {
        vcf_writer_->writeHeader(getMergedHeader());
        sample_indices_ = boost::dynamic_bitset<>(getMergedSamplesNumber());
    }

    ~CombineGVCFs() override {
        ks_free(&tmp_var_ks_);
    }

    void apply(std::vector<std::shared_ptr<VariantContext>> &variants,
        ReferenceContext &reference_context) override;
    
    void finalize() override;
};



#endif  // DPGT_COMBINE_GVCFS_HPP
