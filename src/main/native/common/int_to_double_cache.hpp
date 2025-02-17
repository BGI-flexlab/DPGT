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
#ifndef DPGT_INT_TO_DOUBLE_CACHE_HPP
#define DPGT_INT_TO_DOUBLE_CACHE_HPP

#include <iostream>
#include <cstdlib>


/**
 * @brief abstract base class for int to double cache
 */
class IntToDoubleCache {
private:
    double *cache_ = nullptr;
    int size_ = 0;

public:
    IntToDoubleCache() = default;

    IntToDoubleCache(const IntToDoubleCache &other);
    IntToDoubleCache(IntToDoubleCache &&other) noexcept;

    IntToDoubleCache &operator=(const IntToDoubleCache &other);
    IntToDoubleCache &operator=(IntToDoubleCache &&other) noexcept;

    virtual ~IntToDoubleCache() {
        if (cache_) free(cache_);
    }
    
    virtual int maxSize() = 0;
    virtual double compute(int n) = 0;

    double get(int i);  // not thread safe

    void expandCache(int new_capacity);  // not thread safe

    int size() const {
        return size_;
    }
};



#endif  // DPGT_LOG10_FACTORIAL_CACHE_HPP
