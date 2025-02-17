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
#ifndef DPGT_OVERLAP_DETECTOR_HPP
#define DPGT_OVERLAP_DETECTOR_HPP

#include <map>
#include <vector>
#include "common/simple_interval.hpp"
#include "IITree.h"

/**
 * @brief Detect interval overlap using liheng's implicit interval tree
 * 
 */
class OverlapDetector {
private:
    std::map<int32_t, IITree<int64_t, bool>> interval_map_;
public:
    OverlapDetector() = default;

    void add(const SimpleInterval &interval) {
        if (interval_map_.find(interval.tid) == interval_map_.end()) {
            interval_map_[interval.tid] = IITree<int64_t, bool>();
        }
        interval_map_[interval.tid].add(interval.start, interval.end+1, true);
    }

    void index() {
        for (auto &it: interval_map_) {
            it.second.index();
        }
    }

    std::vector<SimpleInterval> overlap(const SimpleInterval &query) const {
        auto tree_itr = interval_map_.find(query.tid);
        if (tree_itr == interval_map_.end()) return {};
        std::vector<size_t> result_indices;
        tree_itr->second.overlap(query.start, query.end + 1, result_indices);
        std::vector<SimpleInterval> results;
        for (auto i: result_indices) {
            results.emplace_back(query.tid, tree_itr->second.start(i),
                tree_itr->second.end(i) - 1);
        }
        return results;
    }

    bool isOverlap(const SimpleInterval &query) const {
        auto tree_itr = interval_map_.find(query.tid);
        if (tree_itr == interval_map_.end()) return {};
        std::vector<size_t> result_indices;
        tree_itr->second.overlap(query.start, query.end + 1, result_indices);
        return !result_indices.empty();
    }
};




#endif  // DPGT_OVERLAP_DETECTOR_HPP
