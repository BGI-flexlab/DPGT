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
#include "int_to_double_cache.hpp"
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>


IntToDoubleCache::IntToDoubleCache(const IntToDoubleCache &other) {
    this->cache_ = other.cache_;
    this->size_ = other.size_;
}

IntToDoubleCache::IntToDoubleCache(IntToDoubleCache &&other) noexcept {
    this->cache_ = other.cache_;
    this->size_ = other.size_;

    other.cache_ = nullptr;
}

IntToDoubleCache &IntToDoubleCache::operator=(const IntToDoubleCache &other) {
    if (this != &other) {
        // free this
        if (this->cache_) free(this->cache_);

        const size_t bytes = other.size_*sizeof(double);
        this->cache_ = (double *)malloc(bytes);
        this->size_ = other.size_;
        memcpy(this->cache_, other.cache_, bytes);
    }
    return *this;
}

IntToDoubleCache &IntToDoubleCache::operator=(IntToDoubleCache &&other) noexcept
{
    if (this != &other) {
        // free this
        if (this->cache_) free(this->cache_);

        this->cache_ = other.cache_;
        this->size_ = other.size_;

        // set other cache to null
        other.cache_ = nullptr;
    }
    return *this;
}


double IntToDoubleCache::get(int i) {
    if (i < 0) {
        std::cerr << "[IntToDoubleCache] out of range error! "
            "index can not be negative!" << std::endl;
        std::exit(1);
    }
    if (i > size_) {
        if (i > maxSize()) {
            return compute(i);
        }
        int new_capacity = std::max(i + 10, 2 * size_);
        expandCache(new_capacity);
    }
    return cache_[i];
}


void IntToDoubleCache::expandCache(int new_capacity) {
    if (new_capacity < size_) return;
    cache_ = (double *)realloc(cache_, new_capacity*sizeof(double));
    if (cache_ == nullptr) {
        std::cerr << "[IntToDoubleCache] Error! Can not allocate memory!"
            << std::endl;
    }
    for (int i = size_; i < new_capacity; ++i) {
        cache_[i] = compute(i);
    }
}
