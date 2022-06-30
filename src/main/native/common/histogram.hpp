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
