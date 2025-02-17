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
#ifndef DPGT_MULTI_VARIANT_GROUP_ON_START_HPP
#define DPGT_MULTI_VARIANT_GROUP_ON_START_HPP


#include <cstdint>
#include <limits>
#include <vector>
#include <memory>
#include "common/cached_ref_seq.hpp"
#include "common/reference_context.hpp"
#include "htslib/vcf.h"
#include "vcf/multi_vcf_reader.hpp"
#include "vcf/variant_context.hpp"
#include "common/simple_interval.hpp"
#include "common/overlap_detector.hpp"


/**
 * @brief read on multiple vcf files, emit variants which starts at the
 * same position
 */
class MultiVariantGroupOnStart {
protected:
    std::vector<std::shared_ptr<VariantContext>> current_variants_;
    int first_current_var_start_ = 0;
    int last_current_var_start_ = 0;

    MultiVcfReader *reader_ = nullptr;
    CachedRefSeq *cached_ref_seq_ = nullptr;
    std::vector<std::string> merged_samples_;

    // configures

    /**
     * @brief ignore variants that start outside of the target intervals?
     */
    bool ignore_intervals_outside_start_ = false;
    int64_t reference_window_padding_ = 1;


    // target intervals, if intervals_ is nullptr, then assume operate on whole
    // genome regions
    const std::vector<SimpleInterval> *intervals_ = nullptr;

    // for dectect overlap of input interval to target intervals
    OverlapDetector overlap_detector;

    // build overlap detector
    void buildOverlapDetector() {
        if (intervals_ != nullptr) {
            for (auto const &i: *intervals_) {
                overlap_detector.add(i);
            }
            overlap_detector.index();
        }
    }

protected:
    /**
     * @brief test if query locus is within any target interval
     * 
     * @param loc query locus, locus start must equals locus end
     * @return true if query locus is within any target interval
     * @return false query locus is not within any target interval
     */
    bool isWithinInterval(const SimpleInterval &loc) {
        // target interval is not provieded, assume operating on whole genome,
        // so any loc will be within interval
        if (intervals_ == nullptr) return true;
        return overlap_detector.isOverlap(loc);
    }

    ReferenceContext makeSpanningReferenceContext(
        std::vector<std::shared_ptr<VariantContext>> &variants,
        int64_t reference_window_padding) const;
    
    std::vector<std::string> createMergedSamples() {
        bcf_hdr_t *header = getMergedHeader();
        std::vector<std::string> result;
        for (int i = 0; i < bcf_hdr_nsamples(header); ++i) {
            result.push_back(header->samples[i]);
        }
        return result;
    }

public:
    MultiVariantGroupOnStart(MultiVcfReader *reader,
        CachedRefSeq *cached_ref_seq):
        reader_(reader),
        cached_ref_seq_(cached_ref_seq)
    {
        intervals_ = reader_->getIntervals();
        merged_samples_ = createMergedSamples();
        buildOverlapDetector();
    }

    virtual ~MultiVariantGroupOnStart() {}
    
    /**
     * @brief if input variant is start on the same position as current variants
     * add it to current variants, else emit current variants to other function
     * 
     * @param var input variant
     */
    void apply(std::shared_ptr<VariantContext> &&var);

    /**
     * @brief run the pipeline.
     */
    void run();

    virtual void apply(std::vector<std::shared_ptr<VariantContext>> &variants,
        ReferenceContext &reference_context) = 0;
    
    void apply(std::vector<std::shared_ptr<VariantContext>> &variants) {
        ReferenceContext reference_context = makeSpanningReferenceContext(
            variants, reference_window_padding_);
        apply(variants, reference_context);
    }

    const std::vector<SimpleInterval> *getIntervals() const {
        return intervals_;
    }

    // call this on last grouped variants
    virtual void finalize() {
        if (!current_variants_.empty()) apply(current_variants_);
    }

    bcf_hdr_t *getMergedHeader() {
        return reader_->header();
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
};



#endif  // DPGT_MULTI_VARIANT_GROUP_ON_START_HPP
