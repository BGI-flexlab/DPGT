#ifndef SIMPLE_INTERVAL_HPP
#define SIMPLE_INTERVAL_HPP

#include <algorithm>
#include <cstdint>
#include <functional>
#include <vector>
#include <iostream>

#include "common/interval.hpp"
#include "common/math_utils.hpp"
#include "common/sequence_dictionary.hpp"
#include "htslib/sam.h"


/**
 * Class represent a 0-based closed ended genomic interval.
 */
class SimpleInterval: public Interval<int32_t, int64_t> {
public:
    SimpleInterval() = default;
    SimpleInterval(int32_t tid_, int64_t start_, int64_t end_):
        Interval<int32_t, int64_t>(tid_, start_, end_) {}

    bool IsNull() const
    {
        return tid == -1 || start == -1 || end == -1;
    }

    bool operator<(const SimpleInterval &other) const
    {
        return CompareTo(other) < 0;
    }

    bool operator>(const SimpleInterval &other) const
    {
        return CompareTo(other) > 0;
    }

    bool operator==(const SimpleInterval &other) const
    {
        return CompareTo(other) == 0;
    }

    bool operator!=(const SimpleInterval &other) const
    {
        return CompareTo(other) != 0;
    }

    bool operator>=(const SimpleInterval &other) const
    {
        return CompareTo(other) >= 0;
    }

    bool operator<=(const SimpleInterval &other) const
    {
        return CompareTo(other) <= 0;
    }

    int Size() const
    {
        return end - start + 1;
    }

    /**
     * Calculate the overlap size of two interval.
     */
    int Overlap(const SimpleInterval &other) const
    {
        if (tid != other.tid) return 0;

        int ov = std::min(end, other.end) - std::max(start, other.start) + 1;
        return ov > 0 ? ov : 0;
    }

    bool Contain(const SimpleInterval &other) const
    {
        // null region not contain any other region
        if (IsNull()) return false;
        if (tid == other.tid && start <= other.start && end >= other.end) {
            return true;
        }
        return false;
    }

    bool Contain(int32_t other_tid, int64_t other_pos) const
    {
        // null region not contain any other region
        if (IsNull()) return false;
        if (tid == other_tid && start <= other_pos && end >= other_pos) {
            return true;
        }
        return false;
    }

    /**
     * Expand this interval by the length(must be positive), 
     * up bound by target(contig/chromosome) length - 1
     * @param side bit flag indicating expand direction, 1 for left,
     * 2 for right, 3 for both side.
     */
    void Expand(int64_t length, int64_t target_length,
        uint8_t side = 3)
    {
        if ((side & 1) != 0) {
            start -= length;
            // not less than zero
            start = start >= 0 ? start : 0;
        }
        if ((side & 2) != 0) {
            end += length;
            // not exceed target length - 1
            int max_end = target_length - 1;
            end = end <= max_end ? end : max_end;
        }
    }

    /**
     * Expand this interval with respect to other interval.
     * return the merged interval of the two if there are overlap,
     * else return this interval.
     */
    SimpleInterval Expand(const SimpleInterval &other) const
    {
        if (Overlap(other) > 0)
        {
            return (SimpleInterval(
                tid, std::min(start, other.start), std::max(end, other.end)));
        } else {
            return *this;
        }
    }

private:
    int CompareTo(const SimpleInterval &other) const
    {
        if (this == &other)
        {
            return 0;
        } else
        {
            int cmp;
            if ((cmp = Compare(tid, other.tid)) != 0)
            {
                return cmp;
            } else if ((cmp = Compare(start, other.start)) != 0)
            {
                return cmp;
            } else
            {
                return Compare(end, other.end);
            }
        }
    }
};

std::ostream &operator<<(std::ostream &out_stream,
        const SimpleInterval &interval);


/**
 * @brief Merge a given list of SimpleIntervals if overlap. An optional
 * expansion/margin can be used, if nearby intervals are need to be merged.
 * For example:
 *         1----------10          20----------30
 * input         5----------15
 *                             ||
 * output  1----------------15    20----------30
 *
 * @param intervals input chromosome regions
 * @param tid2len callback function for get target length from tid
 * @param exp expand size before merge
 * @return merged intervals
 */
std::vector<SimpleInterval> MergeSimpleInterval(
    const std::vector<SimpleInterval> &intervals,
    std::function<int64_t(int32_t)> tid2len, int64_t exp = 0);


std::vector<SimpleInterval> MergeSimpleInterval(
    const std::vector<SimpleInterval> &intervals,
    const SequenceDictionary &seq_dict, int64_t exp = 0);


std::vector<SimpleInterval> MergeSimpleInterval(
    const std::vector<SimpleInterval> &intervals,
    const sam_hdr_t *head, int64_t exp = 0);


#endif  // SIMPLE_INTERVAL_HPP
