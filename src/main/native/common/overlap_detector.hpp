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
