#include "histogram.hpp"
#include <cstddef>
#include <ios>
#include <limits>
#include <string>
#include <sstream>
#include <iomanip>



const double Histogram::BIN_EPSILON = 0.01;

const double Histogram::NULL_MEDIAN_VALUE = std::numeric_limits<double>::max();


void Histogram::add(double d, int count) {
    int64_t k = getBinnedValue(d);
    if (!isValidBinKey(k)) {
        std::cerr << "[Histogram] Error! " << d << " is outof range."
            << std::endl;
        std::exit(1); 
    }
    auto itr = data_.find(k);
    if (itr == data_.end()) {
        data_[k] = count;
    } else {
        itr->second += count;
    }
}

void Histogram::add(const Histogram &other) {
    if (other.bin_size_ != this->bin_size_) {
        std::cerr << "[Histogram] Error! other histogram bin size is "
            << "not equal with this histogram bin size." << std::endl;
        std::exit(1);
    }
    for (auto &itr: other.data_) {
        auto itr1 = data_.find(itr.first);
        if (itr1 == data_.end()) {
            data_[itr.first] = itr.second;
        } else {
            itr1->second += itr.second;
        }
    }
}

double Histogram::median() const {
    int num_items = 0;
    for (auto &itr: data_) num_items += itr.second;
    bool odd_num_values = num_items % 2 == 0 ? false : true;
    int median_index = (num_items + 1)/2;
    int counter = 0;
    double first_median = NULL_MEDIAN_VALUE;
    for (auto &itr: data_) {
        counter += itr.second;
        if (counter > median_index) {
            if (first_median == NULL_MEDIAN_VALUE) {
                return itr.first * bin_size_;
            } else {
                return (first_median+itr.first)/2.0*bin_size_;
            }
        }
        if (counter == median_index) {
            if (odd_num_values) {
                return itr.first * bin_size_;
            } else {
                first_median = static_cast<double>(itr.first);
            }
        }
    }
    return NULL_MEDIAN_VALUE;
}


int Histogram::get(double d) const {
    int64_t k = getBinnedValue(d);
    if (!isValidBinKey(k)) {
        std::cerr << "[Histogram] Error! " << d << " is outof range."
            << std::endl;
        std::exit(1); 
    }
    auto itr = data_.find(k);
    if (itr == data_.end()) {
        return 0;
    } else {
        return itr->second;
    }
}


std::string Histogram::toString() const {
    if (data_.empty()) return "";
    std::ostringstream ostr;
    size_t n = 0;
    for (auto &itr: data_) {
        ostr << std::fixed << std::setprecision(double_precision_)
            << static_cast<double>(itr.first)*bin_size_ << "," << itr.second;
        ++n;
        if (n != data_.size()) ostr << ",";
    }
    return ostr.str();
}
