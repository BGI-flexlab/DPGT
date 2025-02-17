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
#ifndef DPGT_HISTOGRAM_HPP
#define DPGT_HISTOGRAM_HPP

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <cstdint>
#include <string>


class Histogram {
private:
    std::map<int64_t,int> data_;
    double bin_size_ = 0.1;
    static const double BIN_EPSILON;
    int double_precision_ = 1;

    bool isValidBinKey(int64_t k) const {
        return k>=std::numeric_limits<int>::min() &&
            k<=std::numeric_limits<int>::max();
    }

public:
    static const double NULL_MEDIAN_VALUE;

    Histogram() = default;
    Histogram(double bin_size): bin_size_(bin_size) {
        double_precision_ = std::round(-std::log10(bin_size_));
    }

    void add(double d, int count=1);

    void add(const Histogram &other);

    int get(double d) const;

    double median() const;

    int64_t getBinnedValue(double d) const {
        return static_cast<int64_t>(
            std::floor((d+BIN_EPSILON*bin_size_)/bin_size_));
    }

    std::string toString() const;

    bool empty() const {
        return data_.empty();
    }

};


#endif  // DPGT_HISTOGRAM_HPP
