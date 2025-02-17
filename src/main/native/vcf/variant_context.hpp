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
#ifndef DPGT_VARIANT_CONTEXT_HPP
#define DPGT_VARIANT_CONTEXT_HPP

#include <boost/algorithm/string/join.hpp>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <map>
#include <set>
#include "htslib/vcf.h"
#include "allele.hpp"
#include "genotyper/flat_genotype.hpp"
#include "vcf/vcf_attribute.hpp"
#include "vcf/vcf_id_table.hpp"


enum class VariantContextType:uint8_t {
    NULL_TYPE,
    NO_VARIATION,
    SNP,
    MNP,
    INDEL,
    SYMBOLIC,
    MIXED
};


/**
 * @brief variant context implemented using htslib bcf1_t and bcf_hdr_t
 */
class VariantContext {
public:
    VariantContext() = default;
    VariantContext(bcf1_t *variant, bcf_hdr_t *header, VcfIdTables *id_table);
    VariantContext(const VariantContext &other);
    VariantContext(VariantContext &&other) noexcept;
    VariantContext &operator=(const VariantContext &other);
    VariantContext &operator=(VariantContext &&other) noexcept;
    ~VariantContext();

    int unpack(int which) const {
        return bcf_unpack(variant_, which);
    }

    // accessors

    bool isNull() const {
        return variant_ == nullptr;
    }

    bcf1_t *getVariant() const {
        return variant_;
    }

    bcf_hdr_t *getHeader() const {
        return header_;
    }

    int32_t getContig() const {
        return variant_->rid;
    }

    const char *getContigName() const {
        return bcf_hdr_id2name(header_, variant_->rid);
    }


    /**
     * @brief Get 0-based Start of this variant
     * @return int64_t 
     */
    int64_t getStart() const {
        return variant_->pos;
    }

    /**
     * @brief Get 0-based End of this variant.
     * @return int64_t if "END" in info return info["END"] - 1, else return
     * pos + ref_allele_len - 1
     */
    int64_t getEnd();

    std::string getID() const {
        bcf_unpack(variant_, BCF_UN_STR);
        return variant_->d.id;
    }

    std::vector<std::string> getFilter() const {
        bcf_unpack(variant_, BCF_UN_FLT);
        if (variant_->d.n_flt == 0) return {};
        std::vector<std::string> filters(variant_->d.n_flt);
        for (int i = 0; i < variant_->d.n_flt; ++i) {
            filters[i] = bcf_hdr_int2id(header_, BCF_DT_ID, variant_->d.flt[i]);
        }
        return filters;
    }

    int getNSamples() const {
        return bcf_hdr_nsamples(header_);
    }

    char **getSamples() const {
        return header_->samples;
    }

    const std::vector<int> &getSampleIndices() const {
        return id_tables_->sample_id_table.table;
    }
    
    int getNAllele() const {
        return variant_->n_allele;
    }

    VcfIdTables *getVcfIdTables() const {
        return id_tables_;
    }

    int getMergedNSamples() const {
        return id_tables_->sample_id_table.key_map->size();
    }

    /**
     * @brief Get the Alleles of this variant
     * 
     * @return std::vector<Allele*> vector of allele points
     */
    const std::vector<Allele> &getAlleles();

    /**
     * @brief Get the Reference allele of this variant
     */
    Allele getReference() {
        if (!alleles_.empty()) return alleles_.front();
        if (getNAllele() > 0) alleles_ = createAlleles();
        if (getNAllele() == 0) return Allele();
        return alleles_.front();
    }

    /**
     * Genotype related routines
     */

    const std::vector<FlatGenotype *> &getFlatGenotypes();

    /**
     * @brief Get the INFO fields of this variant
     */
    std::map<std::string, VcfAttributeBase *> &getSharedAttributes();

    VariantContextType getType() {
        if (type_ == VariantContextType::NULL_TYPE) {
            determineType();
        }
        return type_;
    }

    bool isSymbolic() {
        if (type_ == VariantContextType::NULL_TYPE) determineType();
        return type_ == VariantContextType::SYMBOLIC;
    }

    int getMaxPloidy(int default_ploidy=2) {
        int m = 0;
        const std::vector<FlatGenotype *> &genotypes = getFlatGenotypes();
        for (auto &g: genotypes)
        {
            const int p = g->getPloidy();
            if (p>m) m = p;
        }
        if (m == 0) return default_ploidy;
        return m;
    }

    /**
     * @brief get depth(DP) of this variant
     */
    int32_t getDepth() {
        if (DP_ > -1) return DP_;
        DP_ = calculateDepth();
        if (DP_ < 0) DP_ = 0;
        return DP_;
    }

private:
    bcf1_t *variant_ = nullptr;
    bcf_hdr_t *header_ = nullptr;  // shared header
    VcfIdTables *id_tables_ = nullptr;

    int64_t end_ = -1;

    std::vector<Allele> alleles_;

    std::vector<FlatGenotype *> flat_genotypes_;

    std::vector<VcfAttributeGT *> gts_;
    // vcf FORMAT of each sample
    std::vector<std::vector<VcfAttributeBase *>> genotype_attributes_;

    // vcf INFO fields
    std::map<std::string, VcfAttributeBase *> shared_attributes_;

    int32_t DP_ = -1;

    VariantContextType type_ = VariantContextType::NULL_TYPE;

    /**
     * @brief primary genotype(FORMAT of vcf) keys: GT, GQ, DP, AD, PL
     */
    static const std::set<std::string> GENOTYPE_PRIMARY_KEYS;

    void clear();

    /**
     * @brief calculate end position, if there is END field in INFO use it,
     * else calculate it by end = pos + length_of_ref - 1
     */
    int64_t calculateEnd() const;

    /**
     * @brief Create Alleles 
     */
    std::vector<Allele> createAlleles() const;

    /**
     * @brief Create a Genotype Attributes for each sample
     */
    void createGenotypeAttributes();
    
    /**
     * @brief create shared attributes(INFO fields) for this variant
     */
    std::map<std::string, VcfAttributeBase *> createSharedAttributes() const;

    void determineType();

    void determinePolymorphicType();

    VariantContextType typeOfBiallelicVariant(
        const Allele &ref, const Allele &alt);

    std::vector<VcfAttributeBase *> initGenotyeAttributeArray() const {
        return std::vector<VcfAttributeBase *>(
            id_tables_->format_id_table.key_map->size(), nullptr);
    }

    std::vector<std::vector<VcfAttributeBase *>>
    initGenotyeAttributeArrayForSamples() const {
        std::vector<std::vector<VcfAttributeBase *>> result(getNSamples());
        for (int i = 0; i < getNSamples(); ++i) {
            result[i] = initGenotyeAttributeArray();
        }
        return result;
    }

    int32_t calculateDepth();

};



#endif  // DPGT_VARIANT_CONTEXT_HPP
