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
#ifndef DPGT_REFERENCE_CONTEXT_HPP
#define DPGT_REFERENCE_CONTEXT_HPP

#include "common/cached_ref_seq.hpp"
#include "common/simple_interval.hpp"
#include <cstdint>
#include <string>



/**
 * @brief represent local reference bases and reference interval
 */
class ReferenceContext {
private:
    /**
     * @brief backing reference sequence data source
     */
    CachedRefSeq *data_source_ = nullptr;

    /**
     * @brief interval of this reference context
     */
    SimpleInterval interval_;

    /**
     * @brief expanded interval of this reference context
     */
    SimpleInterval window_;

    /**
     * @brief reference sequence of this reference context
     * sequence is from window_
     */
    std::string sequence_;

    int64_t trimToContigStart(int64_t start) const {
        return start < 0 ? 0 : start;
    }

    int64_t trimToContigEnd(int64_t end) const {
        int64_t target_length;
        data_source_->getSequenceDict()->SequenceIndexToLength(interval_.tid,
            target_length);
        return end < target_length ? end : target_length - 1;
    }

public:
    ReferenceContext() = default;
    ReferenceContext(CachedRefSeq *data_source, SimpleInterval interval):
        data_source_(data_source), interval_(std::move(interval)) {}
    

    void setWindow(int64_t leading_bases, int64_t trailing_bases);

    bool isNull() const {
        return interval_.IsNull();
    }

    int32_t getTid() const {
        return interval_.tid;
    }

    int64_t getStart() const {
        return interval_.start;
    }

    int64_t getEnd() const {
        return interval_.end;
    }

    bool hasBackingDataSource() const {
        return data_source_ != nullptr;
    }

    /**
     * @brief Get the local bases of window
     */
    std::string getBases() {
        if (data_source_ == nullptr || window_.IsNull()) {
            return "";
        }

        if (sequence_.empty()) {
            sequence_ = data_source_->getSubsequenceAt(
                window_.tid, window_.start, window_.end);
        }
        return sequence_;
    }

    /**
     * @brief Get the Bases of target interval
     */
    std::string getBases(const SimpleInterval &target_interval) const {
        if (data_source_ == nullptr || target_interval.IsNull()) return "";

        int64_t start = trimToContigStart(target_interval.start);
        int64_t end = trimToContigEnd(target_interval.end);

        return data_source_->getSubsequenceAt(target_interval.tid, start, end);
    }

    /**
     * @brief Get bases from start of interval to end of window
     */
    std::string getForwardBases() {
        std::string bases = getBases();
        int mid = interval_.start - window_.start;  // start index at local sequence
        return bases.substr(mid);
    }

    /**
     * @brief Get interval of this reference context
     */
    const SimpleInterval getInterval() const {
        return interval_;
    }

    /**
     * @brief Get windown of this reference context
     */
    const SimpleInterval getWindow() const {
        return window_;
    }
};



#endif  // DPGT_REFERENCE_CONTEXT_HPP
