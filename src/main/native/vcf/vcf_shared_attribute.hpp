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
#ifndef DPGT_VCF_SHARED_ATTRIBUTE_HPP
#define DPGT_VCF_SHARED_ATTRIBUTE_HPP

#include "common/simple_matrix.hpp"
#include "vcf/vcf_attribute.hpp"
#include "common/histogram.hpp"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <unordered_set>


struct VcfSharedAttributeConstants {
    // separator for allele specific attributes, "|"
    static const std::string ALLELE_SPECIFIC_SEP;

    // keys of attributes that can be reduced(merged)
    static const std::unordered_set<std::string> REDUCIBLE_KEYS;

    // stale keys which will be dropped before combine gvcfs
    static const std::unordered_set<std::string> STALE_KEYS;
};


/**
    * @brief check if attribute key is match AS_*RankSum
    */
bool isASrankSumAttribute(const std::string &key);


/**
 * allele specific rank sum test attributes in INFO field of vcf/bcf
 * For example: AS_RAW_BaseQRankSum, AS_RAW_MQRankSum, AS_RAW_ReadPosRankSum
 */
class ASrankSumAttribute: public VcfAttributeString
{
protected:
    std::vector<Histogram> histogram_data_;

public:

    ASrankSumAttribute(std::string key, int16_t type, int size,
        uint8_t size_type, std::string value):
            VcfAttributeString(std::move(key), type, size,
            size_type, std::move(value)) {}
    
    ASrankSumAttribute(const ASrankSumAttribute &other):
        VcfAttributeString(other)
    {
        histogram_data_ = other.histogram_data_;
    }

    ASrankSumAttribute(ASrankSumAttribute &&other) noexcept:
        VcfAttributeString(other)
    {
        histogram_data_ = std::move(other.histogram_data_);
    }

    ASrankSumAttribute &operator=(const ASrankSumAttribute &other) {
        if (this != &other) {
            VcfAttributeString::operator=(other);
            histogram_data_ = other.histogram_data_;
        }
        return *this;
    }

    ASrankSumAttribute &operator=(ASrankSumAttribute &&other) noexcept {
        if (this != &other) {
            VcfAttributeString::operator=(std::move(other));
            histogram_data_ = std::move(other.histogram_data_);
        }
        return *this;
    }
    
    ~ASrankSumAttribute() override {}

protected:
    /**
     * create histogram from internal string value of this attribute
     */
    void createHistogram();

    std::vector<Histogram> &getHistogram() {
        if (value_.empty()) return histogram_data_;
        if (histogram_data_.empty()) createHistogram();
        return histogram_data_;
    }

    /**
     * set string value for this attribute
     */
    void setValue(std::string value) override {
        value_ = std::move(value);
        histogram_data_.clear();
    }

public:
    /**
    * merge allele specific rank sum attributes by merge histograms
    * @param attributes allele specific rank sum attribute for each variant
    * to be merged
    * @param allele_index_map new to old allele indices map for each variant to be
    * merged, each row is a map of new to old indices
    * @return merged allele specific rank sum attribute
    */
    VcfAttributeBase *merge(
        const std::vector<VcfAttributeBase *> &attributes,
        const SimpleMatrix<int> *allele_index_map_ptr = nullptr) const override;
};


enum class ASAttributeInternalDataType:int8_t {INT, DOUBLE, INVALID};


/**
 * check if string INFO attribute is a common AS attribute and get the internal
 * data type
 */
bool isASAttribute(
    const std::string &key, const std::string &val,
    ASAttributeInternalDataType &type);


/**
 * Class for allele specific attributes of which each sub-field have the same
 * size. For example:
 * AS_RAW_MQ=38|60|40 means that there are 3 alleles, each sub-field(for each
 * allele) have just 1 value.
 * AS_SB_TABLE=12,15|10,2|0,0 means that there are 3 alleles, each sub-field
 * have 2 values.
 * Note that this class not works on AS*RankSum attributes, because sub-fields
 * of them have different size.
 */
template<typename T>
class ASAttribute: public VcfAttributeString {
protected:
    T *data_ = nullptr;
    int data_size_ = 0;
    int sub_field_size_ = 0;

public:
    ASAttribute(std::string key, int16_t type, int size,
        uint8_t size_type, std::string value):
            VcfAttributeString(std::move(key), type, size,
            size_type, std::move(value)) {}
    
    ASAttribute(const ASAttribute &other):
        VcfAttributeString(other)
    {
        this->data_size_ = other.data_size_;
        this->sub_field_size_ = other.sub_field_size_;
        this->data_ = (T *)malloc(this->data_size_*sizeof(T));
        memcpy(this->data_, other.data_, this->data_size_*sizeof(T));
    }

    ASAttribute(ASAttribute &&other) noexcept:
        VcfAttributeString(other)
    {
        this->data_ = other.data_;
        this->data_size_ = other.data_size_;
        this->sub_field_size_ = other.sub_field_size_;
        other.data_ = nullptr;
    }

    ASAttribute &operator=(const ASAttribute &other) {
        if (this != &other) {
            VcfAttributeString::operator=(other);
            free(this->data_);
            this->data_size_ = other.data_size_;
            this->sub_field_size_ = other.sub_field_size_;
            this->data_ = (T *)malloc(this->data_size_*sizeof(T));
            memcpy(this->data_, other.data_, this->data_size_*sizeof(T));
        }
        return *this;
    }

    ASAttribute &operator=(ASAttribute &&other) noexcept {
        if (this != &other) {
            VcfAttributeString::operator=(std::move(other));
            free(this->data_);
            this->data_ = other.data_;
            this->data_size_ = other.data_size_;
            this->sub_field_size_ = other.sub_field_size_;
            other.data_ = nullptr;
        }
        return *this;
    }
    
    ~ASAttribute() override {
        free(this->data_);
    }

protected:
    void createData();

    T *getData(int *data_size, int *sub_field_size) {
        *data_size = data_size_;
        *sub_field_size = sub_field_size_;
        return data_;
    }

    void setValue(std::string value) override {
        value_ = std::move(value);
        free(data_);
        data_size_ = 0;
        sub_field_size_ = 0;
    }

public:
    VcfAttributeBase *merge(
        const std::vector<VcfAttributeBase *> &attributes,
        const SimpleMatrix<int> *allele_index_map_ptr = nullptr) const override;

};


template<typename T>
void ASAttribute<T>::createData() {
    if (value_.empty() || data_ != nullptr) return; 
    std::vector<std::string> fields;
    boost::split(fields, value(),
        boost::is_any_of(
            VcfSharedAttributeConstants::ALLELE_SPECIFIC_SEP));
    
    // get sub field size, assume that every sub field have the same size
    for (auto &f: fields) {
        if (!f.empty() && f != "NaN") {
            std::vector<std::string> sub_fields;
            boost::split(sub_fields, f, boost::is_any_of(","));
            sub_field_size_ = sub_fields.size();
            break;
        }
    }
    
    data_size_ = fields.size()*sub_field_size_;
    data_ = (T *)malloc(data_size_*sizeof(T));

    for (size_t i = 0; i < fields.size(); ++i) {
        int j = i * sub_field_size_;
        auto &f = fields[i];
        if (f.empty() || f == "NaN") {
            for (int k = 0; k < sub_field_size_; ++k) {
                data_[j+k] = static_cast<T>(0); 
            }
            continue;
        }
        std::vector<std::string> sub_fields;
        boost::split(sub_fields, f, boost::is_any_of(","));
        if (static_cast<int>(sub_fields.size()) != sub_field_size_) {
            std::cerr << "[ASAttribute::createData] Error! subfield have "
                << "different size. INFO field: " << key() << "=" << value_
                << std::endl;
            std::exit(1);
        }
        for (int k = 0; k < sub_field_size_; ++k) {
            data_[j+k] = boost::lexical_cast<T>(sub_fields[k]);
        }
    }
}



template<typename T>
VcfAttributeBase * ASAttribute<T>::merge(
    const std::vector<VcfAttributeBase *> &attributes,
    const SimpleMatrix<int> *allele_index_map_ptr) const 
{
    if (allele_index_map_ptr == nullptr) {
        std::cerr << "[ASAttribute::merge] Error! need allele "
            <<"index map to merge." << std::endl;
        std::exit(1);  
    } 
    const SimpleMatrix<int> &allele_index_map = *allele_index_map_ptr;
    if (static_cast<int>(attributes.size()) != allele_index_map.rows()) {
        std::cerr << "[ASAttribute::merge] Error! attribute and "
            << "allele index map not match" << std::endl;
        std::exit(1);
    }

    ASAttribute<T> *first_valid_attribute = nullptr;
    for (auto &a: attributes) {
        if (a != nullptr) {
            first_valid_attribute = dynamic_cast<ASAttribute<T> *>(a);
            break;
        }
    }
    if (first_valid_attribute == nullptr) {
        return nullptr;  // no valid attributes to be merged
    }

    int sub_field_size = -1;
    std::vector<ASAttribute<T> *> as_attributes(
        attributes.size(), nullptr);
    for (size_t i = 0; i < attributes.size(); ++i) {
        if (attributes[i] == nullptr) {
            as_attributes[i] = nullptr;
        } else {
            as_attributes[i] = dynamic_cast<ASAttribute<T> *>(
                attributes[i]);
            as_attributes[i]->createData();
            // check if subfield size are equal
            if (sub_field_size == -1) {
                sub_field_size = as_attributes[i]->sub_field_size_;
            }
            if (as_attributes[i]->sub_field_size_ != sub_field_size) {
                std::cerr << "[ASAttribute::merge] Error! Can not merge "
                    << "attributes which have different subfield size." 
                    << std::endl;
                std::exit(1);
            }
        }
    }

    
    T new_data_size = allele_index_map.cols()*sub_field_size;
    T *new_data = (T *)calloc(new_data_size, sizeof(T));
    
    for (int i = 0; i < allele_index_map.rows(); ++i) {
        // the i-th variant not have this type of attribute
        if (as_attributes[i] == nullptr) continue;
        int old_data_size = 0;
        int old_sub_field_size = 0;
        T *old_data = as_attributes[i]->getData(&old_data_size,
            &old_sub_field_size);
        int p, q, m;
        for (int j = 0; j < allele_index_map.cols(); ++j) {
            p = sub_field_size*j;
            m = sub_field_size*allele_index_map(i, j);
            for (q = 0; q < sub_field_size; ++q) {
                new_data[p+q] += old_data[m+q];
            }
        }
    }

    std::ostringstream new_ostr;
    int last1 = allele_index_map.cols() - 1;
    int last2 = sub_field_size - 1;
    int j;
    for (int i = 0; i < last1; ++i) {
        j = sub_field_size*i;
        for (int k = 0; k < last2; ++k) {
            new_ostr << std::fixed << std::setprecision(2)
                << new_data[j+k] << ",";
        }
        new_ostr << std::fixed << std::setprecision(2) << new_data[j+last2]
            << VcfSharedAttributeConstants::ALLELE_SPECIFIC_SEP;
    }

    j = sub_field_size*last1;
    for (int k = 0; k < last2; ++k) {
        new_ostr << std::fixed << std::setprecision(2) << new_data[j+k] << ",";
    }
    new_ostr << std::fixed << std::setprecision(2) << new_data[j+last2];
    
    free(new_data);

    ASAttribute<T> *new_attribute =
        new ASAttribute<T>(first_valid_attribute->key(),
            first_valid_attribute->type(), first_valid_attribute->size(),
            first_valid_attribute->sizeType(), new_ostr.str());
    return new_attribute;
}


/**
 * RAW_MQandDP attribute, array of two value,
 * first: sum of mq^2, second: total depth
 */
template<typename T>
class VcfSharedAttribute: public VcfAttribute<T> {
public:
    VcfSharedAttribute(std::string key, int16_t type, int size,
        uint8_t size_type, T *value):
        VcfAttribute<T>(std::move(key), type, size, size_type, value) {}
    
    VcfSharedAttribute(const VcfSharedAttribute &other):
        VcfAttribute<T>(other) {}
    
    VcfSharedAttribute(VcfSharedAttribute &&other) noexcept:
        VcfAttribute<T>(std::move(other)) {}
    
    VcfSharedAttribute &operator=(const VcfSharedAttribute &other)
    {
        VcfAttribute<T>::operator=(other);
    }
    
    VcfSharedAttribute &operator=(VcfSharedAttribute &&other)
    {
        VcfAttribute<T>::operator=(std::move(other));
    }

    ~VcfSharedAttribute() override {}

    VcfAttributeBase *merge(
        const std::vector<VcfAttributeBase *> &attributes,
        const SimpleMatrix<int> *allele_index_map_ptr = nullptr) const override;
};


template<typename T>
VcfAttributeBase * VcfSharedAttribute<T>::merge(
    const std::vector<VcfAttributeBase *> &attributes,
    const SimpleMatrix<int> *allele_index_map_ptr) const
{
    int16_t type = -1;
    int size = -1;
    for (auto &itr: attributes) {
        if (itr == nullptr) continue;
        if (type == -1) type = itr->type();
        if (type != itr->type()) {
            std::cerr << "[VcfSharedAttribute::merge] Error! input "
                << "attributes have different types." << std::endl;
            std::exit(1);
        }

        if (itr->type() == BCF_HT_STR) {
            std::cerr << "[VcfSharedAttribute::merge] Error! can not "
                << "merge attributes of string type value." << std::endl;
            std::exit(1);
        }

        if (size == -1) size = itr->size();
        if (size != itr->size()) {
            std::cerr << "[VcfSharedAttribute::merge] Error! input "
                << "attributes have different sizes." << std::endl;
            std::exit(1);
        }
    }

    VcfSharedAttribute<T> *first_valid_attribute = nullptr;
    for (auto &a: attributes) {
        if (a != nullptr) {
            first_valid_attribute =
                dynamic_cast<VcfSharedAttribute<T> *>(a);
            break;
        }
    }
    if (first_valid_attribute == nullptr) {
        return nullptr;  // no valid attributes to be merged
    }

    std::vector<VcfSharedAttribute<T> *> cast_attributes(
        attributes.size(), nullptr);
    for (size_t i = 0; i < attributes.size(); ++i) {
        if (attributes[i] == nullptr) {
            cast_attributes[i] = nullptr;
        } else {
            cast_attributes[i] = dynamic_cast<VcfSharedAttribute<T> *>(
                attributes[i]);
        }
    }

    T *new_value = (T *)calloc(size, sizeof(T));

    for (auto &itr: cast_attributes) {
        if (itr == nullptr) continue;
        for (int i = 0; i < size; ++i) {
            new_value[i] += itr->value()[i];
        }
    }

    VcfSharedAttribute<T> *new_attribute =
        new VcfSharedAttribute<T>(first_valid_attribute->key(),
            first_valid_attribute->type(), first_valid_attribute->size(),
            first_valid_attribute->sizeType(), new_value);
    return new_attribute;
}


#endif  // DPGT_VCF_SHARED_ATTRIBUTE_HPP
