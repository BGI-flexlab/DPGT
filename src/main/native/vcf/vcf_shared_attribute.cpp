#include "vcf/vcf_shared_attribute.hpp"
#include "common/simple_matrix.hpp"
#include "vcf/gatk_vcf_constants.hpp"
#include "vcf/vcf_constants.hpp"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <immintrin.h>


const std::string VcfSharedAttributeConstants::ALLELE_SPECIFIC_SEP = "|";

const std::unordered_set<std::string>
VcfSharedAttributeConstants::REDUCIBLE_KEYS = {
    "AS_RAW_BaseQRankSum", "AS_RAW_MQ", "AS_RAW_MQRankSum", 
    "AS_RAW_ReadPosRankSum", "AS_SB_TABLE", "AS_SOR", "RAW_MQandDP",
    "RAW_MQ"};


const std::unordered_set<std::string> VcfSharedAttributeConstants::STALE_KEYS =
{
    VCFConstants::DEPTH_KEY,  // drop depth key before merge, because depth is calculated separately
    VCFConstants::ALLELE_COUNT_KEY,
    VCFConstants::ALLELE_FREQUENCY_KEY,
    VCFConstants::ALLELE_NUMBER_KEY,
    GATKVCFConstants::MLE_ALLELE_COUNT_KEY,
    GATKVCFConstants::MLE_ALLELE_FREQUENCY_KEY,
    VCFConstants::END_KEY,
    GATKVCFConstants::EVENT_COUNT_IN_HAPLOTYPE_KEY
};



bool isASrankSumAttribute(const std::string &key)
{
    if (boost::starts_with(key, GATKVCFConstants::ALLELE_SPECIFIC_PREFIX) &&
        boost::ends_with(key, "RankSum"))
    {
        return true;
    }
    return false;
}

void ASrankSumAttribute::createHistogram() {
    if (value_.empty() || !histogram_data_.empty()) return; 
    std::vector<std::string> fields;
    boost::split(fields, value(),
        boost::is_any_of(
            VcfSharedAttributeConstants::ALLELE_SPECIFIC_SEP));
    for (auto &f: fields) {
        if (f.empty() || f == "NaN") {
            histogram_data_.push_back(Histogram(0.1));
            continue;
        }
        std::vector<std::string> sub_fields;
        boost::split(sub_fields, f, boost::is_any_of(","));
        // check if sub fields size is even
        if (sub_fields.size() % 2 != 0) {
            std::cerr << "[ASrankSumAttribute::createHistogramData] Error! "
                << "allele specific data size is not even!" << std::endl;
            std::exit(1);
        }
        Histogram h1(0.1);
        for (size_t i = 0, j = 1; j < sub_fields.size(); i+=2, j+=2) {
            h1.add(boost::lexical_cast<double>(sub_fields[i]),
                boost::lexical_cast<int>(sub_fields[j]));
        }
        histogram_data_.push_back(std::move(h1));
    }
}


VcfAttributeBase * ASrankSumAttribute::merge(
        const std::vector<VcfAttributeBase *> &attributes,
        const SimpleMatrix<int> *allele_index_map_ptr) const
{
    if (allele_index_map_ptr == nullptr) {
        std::cerr << "[ASrankSumAttribute::merge] Error! need allele index map"
            << " to merge." << std::endl;
        std::exit(1);  
    } 
    const SimpleMatrix<int> &allele_index_map = *allele_index_map_ptr;
    if (static_cast<int>(attributes.size()) != allele_index_map.rows()) {
        std::cerr << "[ASrankSumAttribute::merge] Error! attribute and allele"
            << " index map not match" << std::endl;
        std::exit(1);
    }

    ASrankSumAttribute *first_valid_attribute = nullptr;
    for (auto &a: attributes) {
        if (a != nullptr) {
            first_valid_attribute = dynamic_cast<ASrankSumAttribute *>(a);
            break;
        }
    } 
    if (first_valid_attribute == nullptr) {
        return nullptr;  // no valid attributes to be merged
    }

    std::vector<ASrankSumAttribute *> as_ranksum_attributes(
        attributes.size(), nullptr);
    for (size_t i = 0; i < attributes.size(); ++i) {
        if (attributes[i] == nullptr) {
            as_ranksum_attributes[i] = nullptr;
        } else {
            as_ranksum_attributes[i] = dynamic_cast<ASrankSumAttribute *>(
                attributes[i]);
        }
    }
    
    std::vector<Histogram> new_hists(
        allele_index_map.cols(), Histogram(0.1));
    
    for (int i = 0; i < allele_index_map.rows(); ++i) {
        // the i-th variant not have this type of attribute
        if (as_ranksum_attributes[i] == nullptr) continue;
        std::vector<Histogram> &old_hists =
            as_ranksum_attributes[i]->getHistogram();
        for (int j = 0; j < allele_index_map.cols(); ++j) {
            new_hists[j].add(old_hists.at(allele_index_map(i, j)));
        }
    }

    std::vector<std::string> new_hists_str;
    new_hists_str.reserve(new_hists.size());
    for (auto const &h: new_hists) new_hists_str.push_back(h.toString());
    ASrankSumAttribute *new_attribute =
        new ASrankSumAttribute(first_valid_attribute->key(),
            first_valid_attribute->type(), first_valid_attribute->size(),
            first_valid_attribute->sizeType(), boost::join(new_hists_str,
            VcfSharedAttributeConstants::ALLELE_SPECIFIC_SEP));
    return new_attribute;
}


bool isASAttribute(
    const std::string &key, const std::string &val,
    ASAttributeInternalDataType &type)
{
    if (boost::starts_with(key, GATKVCFConstants::ALLELE_SPECIFIC_PREFIX))
    {
        // only a roughly check
        if (std::strspn(val.c_str(), "-.0123456789|,") == val.size()) {
            if (val.find('.') != std::string::npos) {
                type = ASAttributeInternalDataType::DOUBLE;
            } else {
                type = ASAttributeInternalDataType::INT;
            }
            return true;
        }
    }
    return false;
}

