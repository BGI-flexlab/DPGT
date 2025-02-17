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
#ifndef DPGT_LOG10_FACTORIAL_CACHE_HPP
#define DPGT_LOG10_FACTORIAL_CACHE_HPP


#include <cmath>
#include <math.h>
#include "int_to_double_cache.hpp"


class log10FactorialCache: public IntToDoubleCache {
private:
    static const int CACHE_SIZE = 10000;

public:
    log10FactorialCache() = default;

    log10FactorialCache(const log10FactorialCache &other):
        IntToDoubleCache(other) {}
    log10FactorialCache(log10FactorialCache &&other) noexcept:
        IntToDoubleCache(std::move(other)) {}

    log10FactorialCache &operator=(const log10FactorialCache &other) noexcept {
        IntToDoubleCache::operator=(other);
        return *this;
    }

    log10FactorialCache &operator=(log10FactorialCache &&other) noexcept {
        IntToDoubleCache::operator=(std::move(other));
        return *this;
    }

    ~log10FactorialCache() override {}
    
    int maxSize() override {
        return CACHE_SIZE;
    }

    double compute(int n) override {
        static const double LOG10_E = std::log10(std::exp(1));
        return std::lgamma(n+1) * LOG10_E;
    }

};



#endif  // DPGT_LOG10_FACTORIAL_CACHE_HPP
