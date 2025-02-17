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
#include "simple_interval.hpp"
#include "htslib/sam.h"
#include <cstdint>


std::vector<SimpleInterval> MergeSimpleInterval(
    const std::vector<SimpleInterval> &intervals,
    std::function<int64_t(int32_t)> tid2len, int64_t exp)
{
    std::vector<SimpleInterval> merged_intervals;
    std::vector<const SimpleInterval *> intervals_ptr;
    for (auto const &it: intervals)
    {
        intervals_ptr.push_back(&it);
    }

    std::sort( intervals_ptr.begin(), intervals_ptr.end(), 
        [] (const SimpleInterval *a, const SimpleInterval *b)
        { return *a < *b; } );
    
    SimpleInterval pre;
    for (auto const &r: intervals_ptr)
    {
        SimpleInterval cur = *r;
        cur.Expand(exp, tid2len(r->tid));
        if (pre.IsNull())
        {
            pre = cur;
            continue;
        }
        if (cur.Overlap(pre) > 0)
        {
            pre = pre.Expand(cur);
        } else
        {
            merged_intervals.push_back(pre);
            pre = cur;
        }
    }
    if (!pre.IsNull()) merged_intervals.push_back(pre);
    return merged_intervals;
}


std::vector<SimpleInterval> MergeSimpleInterval(
    const std::vector<SimpleInterval> &intervals,
    const SequenceDictionary &seq_dict, int64_t exp)
{
    auto tid2len = [&seq_dict](int32_t tid) -> int64_t {
        int64_t len;
        seq_dict.SequenceIndexToLength(tid, len);
        return len;
    };

    return MergeSimpleInterval(intervals, tid2len, exp);
}


std::vector<SimpleInterval> MergeSimpleInterval(
    const std::vector<SimpleInterval> &intervals,
    const sam_hdr_t *head, int64_t exp)
{
    auto tid2len = [head](int32_t tid) -> int64_t {
        return sam_hdr_tid2len(head, tid);
    }; 

    return MergeSimpleInterval(intervals, tid2len, exp);
}


std::ostream &operator<<(std::ostream &out_stream,
        const SimpleInterval &interval)
{
    out_stream << interval.tid << ":" << interval.start << "-" << interval.end;
    return out_stream;
}

