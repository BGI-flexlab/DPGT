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
#ifndef DPGT_VCF_ATTRIBUTE_HPP
#define DPGT_VCF_ATTRIBUTE_HPP

#include <array>
#include <cstdint>
#include <cstdio>
#include <ios>
#include <string>
#include <map>
#include <sstream>
#include <type_traits>
#include <vector>
#include <iomanip>
#include <cmath>
#include "htslib/vcf.h"
#include "vcf/vcf_constants.hpp"
#include "common/simple_matrix.hpp"

/**
 * @brief represent attribute of vcf FORMAT or info
 */
class VcfAttributeBase {
protected:
    std::string key_;

    /**
     * @brief type of the attribute value as defined in htslib vcf.h
     * 
     * #define BCF_HT_FLAG 0 // header type
     * #define BCF_HT_INT  1
     * #define BCF_HT_REAL 2
     * #define BCF_HT_STR  3
     * #define BCF_HT_LONG (BCF_HT_INT | 0x100) // BCF_HT_INT, but for int64_t values; VCF only!
     * 
     */
    int16_t type_;

    // string explain of types
    static const std::map<int, std::string> TYPE_STRINGS;

    // number of elements for the value
    int size_;

    /**
     * @brief type of the number of elements as defined in htslib vcf.h
     * 
     * #define BCF_VL_FIXED 0
     * #define BCF_VL_VAR   1 // variable length
     * #define BCF_VL_A     2
     * #define BCF_VL_G     3
     * #define BCF_VL_R     4
     */
    uint8_t size_type_;

    static const std::array<std::string, 5> SIZE_TYPE_STRINGS;

public:
    VcfAttributeBase(std::string key, int16_t type, int size,
        uint8_t size_type):
        key_(std::move(key)), type_(type), size_(size), size_type_(size_type) {}
    
    VcfAttributeBase(const VcfAttributeBase &other) = default;
    VcfAttributeBase(VcfAttributeBase &&other) noexcept {
        key_ = other.key_;
        type_ = other.type_;
        size_ = other.size_;
        size_type_ = other.size_;
    }

    VcfAttributeBase &operator=(const VcfAttributeBase &other) {
        if (this != &other) {
            key_ = other.key_;
            type_ = other.type_;
            size_ = other.size_;
            size_type_ = other.size_;
        }
        return *this;
    }

    VcfAttributeBase &operator=(VcfAttributeBase &&other) noexcept {
        if (this != &other) {
            key_ = std::move(other.key_);
            type_ = other.type_;
            size_ = other.size_;
            size_type_ = other.size_;
        }
        return *this;
    }

    virtual ~VcfAttributeBase() {}

    std::string key() const {
        return key_;
    }

    int type() const {
        return type_;
    }

    const std::string &typeString() const {
        return TYPE_STRINGS.at(type_);
    }

    int size() const {
        return size_;
    }

    uint8_t sizeType() const {
        return size_type_;
    }

    const std::string &sizeTypeString() const {
        return SIZE_TYPE_STRINGS[size_type_];
    }

    virtual void getValueStr(kstring_t *s) = 0;

    virtual VcfAttributeBase *merge(
        const std::vector<VcfAttributeBase *> &attributes,
        const SimpleMatrix<int> *allele_index_map_ptr = nullptr) const
    {
        // do nothing
        return nullptr;
    }
};


template <typename T>
class VcfAttribute: public VcfAttributeBase {
protected:
    T *value_;

    void formatVCFDouble(kstring_t *s, double d);
public:
    VcfAttribute(std::string key, int16_t type, int size,
        uint8_t size_type, T *value):
        VcfAttributeBase(std::move(key), type, size, size_type),
        value_(value) {}
    
    VcfAttribute(const VcfAttribute &other): VcfAttributeBase(other) {
        this->value_ = (T *)malloc(this->size_*sizeof(T));
        memcpy(this->value_, other.value_, this->size_*sizeof(T));
    }

    VcfAttribute(VcfAttribute &&other) noexcept:
        VcfAttributeBase(std::move(other))
    {
        this->value_ = other.value_;
        other.value_ = nullptr;
    }

    VcfAttribute &operator=(const VcfAttribute &other) {
        if (this != &other) {
            clear();
            VcfAttributeBase::operator=(other);
            this->value_ = (T *)malloc(this->size_*sizeof(T));
            memcpy(this->value_, other.value_, this->size_*sizeof(T));
        }
        return *this;
    }

    VcfAttribute &operator=(VcfAttribute &&other) noexcept {
        if (this != &other) {
            clear();
            VcfAttributeBase::operator=(std::move(other));
            this->value_ = other.value_;
            other.value_ = nullptr;
        }
        return *this;
    }

    ~VcfAttribute() override {
        clear();
    }

    T &operator[](int i) {
        return value_[i];
    }

    const T &operator[](int i) const {
        return value_[i];
    }

    T *value() const {
        return value_;
    }

    void clear() {
        if (value_) free(value_);
    }

    void getValueStr(kstring_t *s) override {
        return getValueStrImpl(s);
    }

    template<typename Q = T>
    typename std::enable_if<std::is_floating_point<Q>::value,
        void>::type
    getValueStrImpl(kstring_t *s) {
        if (size_ == 0) return;
        const int last = size_ - 1;
        for (int i = 0; i < last; ++i) {
            formatVCFDouble(s, value_[i]);
            kputc(',', s);
        }
        formatVCFDouble(s, value_[last]);
    }

    template<typename Q = T>
    typename std::enable_if<!std::is_floating_point<Q>::value,
        void>::type
    getValueStrImpl(kstring_t *s) {
        if (size_ == 0) return;
        const int last = size_ - 1;
        for (int i = 0; i < last; ++i) {
            kputl(value_[i], s);
            kputc(',', s);
        }
        kputl(value_[last], s);
    }

// not use the c++17 feature 'if constexpr'
#if 0
    std::string getValueStrImpl() const {
        if (size_ == 0) return "";
        std::ostringstream out_str;
        const int last = size_ - 1;
        if constexpr(std::is_floating_point<T>::value) {
            for (int i = 0; i < last; ++i) {
                out_str << formatVCFDouble(value_[i]) << ",";
            }
        } else {
            for (int i = 0; i < last; ++i) {
                out_str << value_[i] << ",";
            }
        }
        if constexpr(std::is_floating_point<T>::value) {
            out_str << formatVCFDouble(value_[last]);
        } else {
            out_str << value_[last];
        }
        return out_str.str();
    }
#endif
};

template<typename T>
void VcfAttribute<T>::formatVCFDouble(kstring_t *s, double d) {
    if (d < 1.0) {
        if (d < 0.01) {
            if (fabs(d) < 1.0e-20) {
                kputs("0.00", s);
                return;
            }
            ksprintf(s, "%.3e", d);
        } else {
            ksprintf(s, "%.3f", d);
        }
    } else {
        ksprintf(s, "%.2f", d);
    }
}


class VcfAttributeString: public VcfAttributeBase {
protected:
    std::string value_;

public:
    VcfAttributeString(std::string key, int16_t type, int size,
        uint8_t size_type, std::string value):
        VcfAttributeBase(std::move(key), type, size, size_type),
        value_(std::move(value)) {}
    
    VcfAttributeString(const VcfAttributeString &other): VcfAttributeBase(other)
    {
        value_ = other.value_;
    }

    VcfAttributeString(VcfAttributeString &&other) noexcept:
        VcfAttributeBase(std::move(other))
    {
        this->value_ = std::move(other.value_);
    }

    VcfAttributeString &operator=(const VcfAttributeString &other) {
        if (this != &other) {
            VcfAttributeBase::operator=(other);
            value_ = other.value_;
        }
        return *this;
    }

    VcfAttributeString &operator=(VcfAttributeString &&other) noexcept {
        if (this != &other) {
            VcfAttributeBase::operator=(std::move(other));
            this->value_ = std::move(other.value_);
        }
        return *this;
    }

    ~VcfAttributeString() override {
    }

    const std::string &value() const {
        return value_;
    }

    virtual void setValue(std::string new_value) {
        value_ = std::move(new_value);
    }

    void getValueStr(kstring_t *s) override {
        if (size_ == 0) return;
        kputs(value_.c_str(), s);
    }
};


class VcfAttributeGT: public VcfAttribute<int32_t> {
protected:
    bool is_phased_;
public:
    VcfAttributeGT(std::string key, int16_t type, int size, uint8_t size_type,
        int32_t *value):
        VcfAttribute<int32_t>(std::move(key), type, size, size_type, value)
    {
        for (int i = 0; i < size_; ++i) {
            if (value_[i] == bcf_int32_missing) {
                continue;
            } else {
                const int32_t a = bcf_gt_allele(value_[i]);
                is_phased_ = bcf_gt_is_phased(value_[i]);
                value_[i] = a;
            }
        }
    }

    VcfAttributeGT(const VcfAttributeGT &other): VcfAttribute<int32_t>(other)
    {
        this->is_phased_ = other.is_phased_;
    }

    VcfAttributeGT(VcfAttributeGT &&other) noexcept:
        VcfAttribute<int32_t>(std::move(other))
    {
        this->is_phased_ = other.is_phased_;
    }

    VcfAttributeGT &operator=(const VcfAttributeGT &other) {
        if (this != &other) {
            VcfAttribute<int32_t>::operator=(other);
            this->is_phased_ = other.is_phased_;
        }
        return *this;
    }

    VcfAttributeGT &operator=(VcfAttributeGT &&other) noexcept {
        if (this != &other) {
            VcfAttribute<int32_t>::operator=(std::move(other));
            this->is_phased_ = other.is_phased_;
        }
        return *this;
    }

    ~VcfAttributeGT() override {}

    bool isPhased() const {
        return is_phased_;
    }

    // set GT to NO_CALL/MISSING by assign bcf_int32_missing to values

    /**
     * @brief Set GT to NO_CALL/MISSING by set values to bcf_int32_missing
     * 
     */
    void setToNoCall() {
        for (int i = 0; i < size_; ++i) {
            value_[i] = bcf_int32_missing;
        }
    }

    void getValueStr(kstring_t *s) override {
        if (size_ == 0) return;
        const int last = size_ - 1;
        for (int i = 0; i < last; ++i) {
            if (value_[i] < 0) {
                kputc('.', s);
            } else {
                kputl(value_[i], s);
            }
            if (is_phased_) {
                kputs(VCFConstants::PHASED.c_str(), s);
            } else {
                kputs(VCFConstants::UNPHASED.c_str(), s);
            }
        }
        if (value_[last] < 0) {
            kputc('.', s);
        } else {
            kputl(value_[last], s);
        }
    }

};


#endif  // DPGT_VCF_ATTRIBUTE_HPP
