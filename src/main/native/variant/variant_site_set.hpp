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
#ifndef DPGT_VARIANT_SITE_SET_HPP
#define DPGT_VARIANT_SITE_SET_HPP

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include "boost/dynamic_bitset/dynamic_bitset.hpp"

/**
 * Represent variant sites of a genome region using dynamic bitset
 */
class VariantSiteSet {
private:
    boost::dynamic_bitset<uint8_t> data_;
    int64_t start_;  // 0-based start
    int64_t end_;    // 0-based end
    int64_t size_;

public:
    VariantSiteSet() = default;

    VariantSiteSet(int64_t start, int64_t end): start_(start), end_(end) {
        size_ = end - start + 1;
        data_ = boost::dynamic_bitset<uint8_t>(size_);
    }

    VariantSiteSet &set(int64_t pos) {
        validatePosition(pos);
        data_.set(pos - start_);
        return *this;
    }

    VariantSiteSet &reset(int64_t pos) {
        validatePosition(pos);
        data_.reset(pos - start_);
        return *this;
    }

    bool test(int64_t pos) const {
        validatePosition(pos);
        return data_.test(pos - start_);
    }

    bool testIndex(int64_t i) const {
        return data_.test(i);
    }

    int64_t get(int64_t i) const {
        return i + start_;
    }

    int64_t size() const {
        return size_;
    }


    boost::dynamic_bitset<uint8_t> &data() {
        return data_;
    }

    const boost::dynamic_bitset<uint8_t> &data() const {
        return data_;
    }

    int64_t firstSite() const {
        int i = 0;
        for (; i < size_; ++i) {
            if (testIndex(i)) {
                return i + start_;
            }
        }
        return -1;  // first site not exist
    }

    int64_t nextSite(int64_t &site) const {
        validatePosition(site);
        int i = site - start_ + 1;
        for (; i < size_; ++i) {
            if (testIndex(i)) {
                return i + start_;
            }
        }
        return -1;  // next site not exist
    }

    std::vector<uint8_t> toUint8Vector() const
    {
        std::vector<uint8_t> bytes;
        bytes.reserve(data_.num_blocks());
        boost::to_block_range(data_, std::back_inserter(bytes));
        return bytes;
    }

    void fromUnit8Vector(const std::vector<uint8_t> &bytes) {
        boost::from_block_range(bytes.begin(), bytes.end(), data_);
    }

private:
    void validatePosition(int64_t pos) const {
        if (pos < start_ || pos > end_) {
            std::cerr << "[VariantSiteSet] Error! Position" << pos
                << " out of range[" << start_ << "," << end_ << "].";
            std::exit(1);
        }
    }
};

#endif  // DPGT_VARIANT_SITE_SET_HPP
