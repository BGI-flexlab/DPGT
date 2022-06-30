#ifndef DPGT_VARIANT_BUILDER_HPP
#define DPGT_VARIANT_BUILDER_HPP

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
#include "htslib/vcf.h"
#include "vcf/allele.hpp"
#include "vcf/vcf_attribute.hpp"
#include "genotyper/flat_genotype.hpp"
#include "vcf/vcf_constants.hpp"


/**
 * @brief build vcf record
 */
class VariantBuilder {
private:
    static const float FLOAT_MISSING;
    bcf_hdr_t *header_;
    VcfKeyMaps *key_maps_;
    int32_t tid_;
    int64_t start_;
    std::string id_ = VCFConstants::MISSING_VALUE_v4;
    int64_t end_;
    std::vector<Allele> alleles_;
    float qual_;
    std::string filter_;
    std::map<std::string, VcfAttributeBase *> *attributes_;
    std::vector<FlatGenotype *> *genotypes_;

    int getMaxPloidy(int default_ploidy = 2);

public:
    VariantBuilder(bcf_hdr_t *header, VcfKeyMaps *key_maps,
        int32_t tid, int64_t start, int64_t end, std::vector<Allele> alleles);

    VariantBuilder &setID(std::string id) {
        id_ = std::move(id);
        return *this;
    }

    VariantBuilder &setFilter(std::string filter) {
        filter_ = std::move(filter);
        return *this;
    }
    
    VariantBuilder &setAttributes(
        std::map<std::string, VcfAttributeBase *> *attributes);

    VariantBuilder &setGenotypes(
        std::vector<FlatGenotype *> *genotypes);

    std::map<std::string, VcfAttributeBase *> &getAttributes() {
        return *attributes_;
    }

    std::vector<FlatGenotype *> &getGenotypes() {
        return *genotypes_;
    }

    /**
     * @brief make variant record string
     */
    void makeString(kstring_t *out_ks);

    int32_t tid() const { return tid_; }
    int64_t start() const { return start_; }
    int64_t end() const { return end_; }
};


struct VariantStringWithPos {
    std::string var_str;
    int32_t rid;
    int64_t start;
    int64_t end;
};


#endif  // DPGT_VARIANT_BUILDER_HPP
