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
#include "variant_context.hpp"
#include "common/math_utils.hpp"
#include "genotyper/genotype.hpp"
#include "htslib/vcf.h"
#include "vcf/allele.hpp"
#include "vcf/vcf_attribute.hpp"
#include "vcf/vcf_shared_attribute.hpp"
#include "gatk_vcf_constants.hpp"
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <vector>


const std::set<std::string> VariantContext::GENOTYPE_PRIMARY_KEYS = {
    "GT", "GQ" ,"DP", "AD", "PL"};


VariantContext::VariantContext(bcf1_t *variant, bcf_hdr_t *header,
    VcfIdTables *id_tables) {
    variant_ = bcf_dup(variant);
    header_ = header;
    id_tables_ = id_tables;
}

VariantContext::VariantContext(const VariantContext &other) {
    this->variant_ = bcf_dup(other.variant_);
    this->header_ = other.header_;
    this->id_tables_ = other.id_tables_;

    this->end_ = other.end_;
    this->alleles_ = other.alleles_;
}

VariantContext::VariantContext(VariantContext &&other) noexcept {
    this->variant_ = other.variant_;
    this->header_ = other.header_;
    this->id_tables_ = other.id_tables_;

    this->end_ = other.end_;
    this->alleles_ = std::move(other.alleles_);

    other.variant_ = nullptr;
    other.header_ = nullptr;
    
    other.alleles_.clear();
}

VariantContext &VariantContext::operator=(const VariantContext &other) {
    if (&other != this) {
        clear();
        
        this->variant_ = bcf_dup(other.variant_);
        this->header_ = other.header_;
        this->id_tables_ = other.id_tables_;
        
        this->end_ = other.end_;
        this->alleles_ = other.alleles_;
    }
    return *this;
}

VariantContext &VariantContext::operator=(VariantContext &&other) noexcept {
    if (&other != this) {
        clear();

        this->variant_ = other.variant_;
        this->header_ = other.header_;
        this->id_tables_ = other.id_tables_;

        this->end_ = other.end_;
        this->alleles_ = std::move(other.alleles_);

        other.variant_ = nullptr;
        other.header_ = nullptr;
        
        other.alleles_.clear();
    }
    return *this;
}

VariantContext::~VariantContext() {
    clear();
}

int64_t VariantContext::getEnd() {
    if (end_ > -1) return end_;
    end_ = calculateEnd();
    return end_;
}

const std::vector<Allele> &VariantContext::getAlleles() {
    if (!alleles_.empty()) return alleles_;
    if (getNAllele() > 0) alleles_ = createAlleles();
    return alleles_;
}

const std::vector<FlatGenotype *> &VariantContext::getFlatGenotypes()
{   
    if (!flat_genotypes_.empty()) return flat_genotypes_;
    if (gts_.empty() && genotype_attributes_.empty()) {
        createGenotypeAttributes();
    }

    flat_genotypes_ = std::vector<FlatGenotype *>(
        genotype_attributes_.size(), nullptr);
    // loop through each sample
    for (size_t i = 0; i < genotype_attributes_.size(); ++i) {
        flat_genotypes_[i] = new FlatGenotype(
            id_tables_->sample_id_table.table[i], gts_[i], genotype_attributes_[i],
            dynamic_cast<VcfFormatKeyMap *>(id_tables_->format_id_table.key_map));
    }

    return flat_genotypes_;
}

std::map<std::string, VcfAttributeBase *> &
VariantContext::getSharedAttributes()
{
    if (!shared_attributes_.empty()) return shared_attributes_;
    shared_attributes_ = createSharedAttributes();
    return shared_attributes_;
}

int64_t VariantContext::calculateEnd() const {
    int end;
    int n = 0;
    int64_t *dst = NULL;
    int ret = bcf_get_info_int64(header_, variant_, "END", &dst, &n);
    if (ret <= 0) {
        // can not find "END" in info
        // end = pos + ref_allele_len - 1
        end = variant_->pos + variant_->rlen - 1;
        return end;
    } else {
        end = *dst - 1;
        if (end < 0) end = 0;
        free(dst);
        return end;
    }
}

std::vector<Allele> VariantContext::createAlleles() const {
    if (!(variant_->unpacked&BCF_UN_STR)) bcf_unpack(variant_, BCF_UN_STR);
    std::vector<Allele> alleles(getNAllele());
    if (getNAllele() > 0) {
        alleles[0] = Allele::create(variant_->d.allele[0], true);
    }
    for (int i = 1; i < getNAllele(); ++i) {
        alleles[i] = Allele::create(variant_->d.allele[i], false);
    }
    return alleles;
}


int32_t VariantContext::calculateDepth() {
    const std::map<std::string, VcfAttributeBase *> &shared_attributes =
        getSharedAttributes();
    auto itr = shared_attributes.find(VCFConstants::DEPTH_KEY);

    // get DP from INFO
    if (itr != shared_attributes.end()) {
        if (itr->second->type() != BCF_HT_INT) {
            std::cerr << "[ReferenceConfidentVariantMerger::calculateVCDepth] "
                << "Error! DP field of vcf INFO must have an Integer value."
                << std::endl;
            return 0;
        }
        VcfAttribute<int32_t> *dp = dynamic_cast<VcfAttribute<int32_t> *>(
            itr->second);
        if (dp->value()[0] > 0) {
            return dp->value()[0];
        }
    }

    // can not get MIN_DP/DP from INFO, so calculate it from genotypes of each sample
    int dp = 0;
    std::vector<FlatGenotype *> genotypes = getFlatGenotypes();
    for (auto &g: genotypes) {
        VcfAttribute<int32_t> *min_dp = g->getMIN_DP();
        if (min_dp) {
            dp += min_dp->value()[0];
        } else if (g->hasDP() && g->getDP()->value()[0] > 0) {
            dp += g->getDP()->value()[0];
        } else if (g->hasAD()) {
            for (int i = 0; i < g->getAD()->size(); ++i) {
                if (g->getAD()->value()[i] > 0) dp += g->getAD()->value()[i];
            }
        }
    }

    return dp;
}


void VariantContext::createGenotypeAttributes()
{
    if (!(variant_->unpacked&BCF_UN_FMT)) bcf_unpack(variant_, BCF_UN_FMT);
    const int nsamples = bcf_hdr_nsamples(header_);

    // init GTs
    gts_ = std::vector<VcfAttributeGT *>(getNSamples(), nullptr);    
    genotype_attributes_ = initGenotyeAttributeArrayForSamples();

    // loop through all genotype fields
    for (int i = 0; i < variant_->n_fmt; ++i) {
        // new format id is index in VcfFmtKeyMap
        int new_id = id_tables_->format_id_table.table[variant_->d.fmt[i].id];
        const char *key =
            header_->id[BCF_DT_ID][variant_->d.fmt[i].id].key;
        const int val_type = variant_->d.fmt[i].type;
        switch (val_type) {
            case BCF_BT_INT8:
            case BCF_BT_INT16:
            case BCF_BT_INT32:
            case BCF_BT_INT64:
            {
                // dst , ndst
                int ndst = 0;
                int32_t *dst = nullptr;
                int n = bcf_get_format_int32(
                    header_, variant_, key, &dst, &ndst);
                if (n > 0) {
                    n /= nsamples;
                    for (int j = 0; j < nsamples; ++j) {
                        int32_t *ptr = dst + j*n;
                        if (strcmp(key, "GT") != 0 && *ptr == bcf_int32_missing) {
                            continue;  // skip missing value
                        }
                        int k;
                        for (k = 0; k < n; ++k) {
                            if (ptr[k] == bcf_int32_vector_end) break;
                        }
                        const size_t s = k*sizeof(int32_t);
                        int32_t *tmp = (int32_t *)malloc(s);
                        memcpy(tmp, ptr, s);
                        const uint8_t l = bcf_hdr_id2length(
                            header_, BCF_HL_FMT, variant_->d.fmt[i].id);
                        if (strcmp(key, "GT") == 0) {
                            // GT coded phase into int32_t, so decoded it 
                            // here
                            gts_[j] = new VcfAttributeGT(
                                key, BCF_HT_INT, k, l, tmp);
                        } else {
                            genotype_attributes_[j][new_id] =
                                new VcfAttribute<int32_t>(
                                    key, BCF_HT_INT, k, l, tmp);
                        }
                    }
                }
                free(dst);
                break;
            }
            case BCF_BT_FLOAT:
            {
                // dst , ndst
                int ndst = 0;
                float *dst = nullptr;
                int n = bcf_get_format_float(
                    header_, variant_, key, &dst, &ndst);
                if (n > 0) {
                    n /= nsamples;
                    for (int j = 0; j < nsamples; ++j) {
                        float *ptr = dst + j*n;
                        if (*ptr == bcf_float_missing) {
                            continue;  // skip missing value
                        }
                        int k;
                        for (k = 0; k < n; ++k) {
                            if (ptr[k] == bcf_float_vector_end) break;
                        }
                        const size_t s = k*sizeof(float);
                        float *tmp = (float *)malloc(s);
                        memcpy(tmp, ptr, s);
                        const uint8_t l = bcf_hdr_id2length(
                            header_, BCF_HL_FMT, variant_->d.fmt[i].id);
                        genotype_attributes_[j][new_id] =
                            new VcfAttribute<float>(
                                key, BCF_HT_REAL, k, l, tmp);
                    }
                }
                free(dst);
                break;
            }
            case BCF_BT_CHAR:
            {
                // dst , ndst
                int ndst = 0;
                char **dst = nullptr;
                int n = bcf_get_format_string(
                    header_, variant_, key, &dst, &ndst);
                if (n > 0) {
                    for (int j = 0; j < nsamples; ++j) {
                        const uint8_t l = bcf_hdr_id2length(
                            header_, BCF_HL_FMT, variant_->d.fmt[i].id);
                        genotype_attributes_[j][new_id] =
                            new VcfAttributeString(
                                key, BCF_HT_STR, 1, l, dst[j]);
                    }
                }
                free(dst[0]); free(dst);
                break;
            }
            default:
            {
                std::cerr << "Error! Invalid format value type for key="
                    << key << " value type=" << val_type << std::endl;
                std::exit(1);
            }
        }
    }
}

std::map<std::string, VcfAttributeBase *>
VariantContext::createSharedAttributes() const
{
    if(!(variant_->unpacked&BCF_UN_SHR)) bcf_unpack(variant_, BCF_UN_SHR);

    std::map<std::string, VcfAttributeBase *> result;

    // loop through all info fields
    for (int i = 0; i < variant_->n_info; ++i) {
        const char *key =
            header_->id[BCF_DT_ID][variant_->d.info[i].key].key;
        const int val_type = variant_->d.info[i].type;
        switch (val_type) {
            case BCF_BT_INT8:
            case BCF_BT_INT16:
            case BCF_BT_INT32:
            case BCF_BT_INT64:
            {
                int ndst = 0;
                int32_t *dst = nullptr;
                int n = bcf_get_info_int32(header_, variant_, key, &dst, &ndst);
                if (n >= 0) {
                    const uint8_t l = bcf_hdr_id2length(header_, BCF_HL_INFO,
                        variant_->d.info[i].key);
                    result[key] = new VcfSharedAttribute<int32_t>(
                        key, BCF_HT_INT, ndst, l, dst);
                }
                break;
            }
            case BCF_BT_FLOAT:
            {
                int ndst = 0;
                float *dst = nullptr;
                int n = bcf_get_info_float(header_, variant_, key, &dst, &ndst);
                if (n >= 0) {
                    const uint8_t l = bcf_hdr_id2length(header_, BCF_HL_INFO,
                        variant_->d.info[i].key);
                    result[key] = new VcfSharedAttribute<float>(
                        key, BCF_HT_REAL, ndst, l, dst);
                }
                break;
            }
            case BCF_BT_CHAR:
            {
                // note that, only support attribute value number=1, this is OK
                // for GATK vcf, but need to support other value numbers in the
                // future
                int ndst = 0;
                char *dst = nullptr;
                int n = bcf_get_info_string(
                    header_, variant_, key, &dst, &ndst);
                if (n >= 0) {
                    const uint8_t l = bcf_hdr_id2length(header_, BCF_HL_INFO,
                        variant_->d.info[i].key);
                    ASAttributeInternalDataType as_attribute_internal_type =
                        ASAttributeInternalDataType::INVALID;
                    if (isASrankSumAttribute(key)) {
                        result[key] = new ASrankSumAttribute(
                            key, BCF_HT_STR, 1, l, dst);
                    } else if (
                        isASAttribute(key, dst, as_attribute_internal_type))
                    {
                        if (as_attribute_internal_type ==
                            ASAttributeInternalDataType::DOUBLE)
                        {
                            result[key] = new ASAttribute<double>(
                                key, BCF_HT_STR, 1, l, dst);
                        } else {  // INT
                            result[key] = new ASAttribute<int>(
                                key, BCF_HT_STR, 1, l, dst);
                        }
                    } else {
                        result[key] = new VcfAttributeString(
                            key, BCF_HT_STR, 1, l, dst);
                    }
                }
                free(dst);
                break;
            }
            default:
            {
                break;
                // std::cerr << "Error! Invalid format value type for key="
                //     << key << " value type=" << val_type << std::endl;
                // std::exit(1);
            }
        }
    }
    return result;
}


void VariantContext::clear() {
    if (variant_) {
        bcf_destroy1(variant_);
    }
    if (!flat_genotypes_.empty()) {
        for (auto g: flat_genotypes_) delete g;
        flat_genotypes_.clear();
    }
    if (!gts_.empty()) {
        for (auto &it: gts_) {
            delete it;
        }
        gts_.clear();
    }
    if (!genotype_attributes_.empty()) {
        for (auto &it: genotype_attributes_) {
            for (auto &it1: it) {
                if (it1) delete it1;
            }
        }
        genotype_attributes_.clear();
    }
    if (!shared_attributes_.empty()) {
        for (auto &it: shared_attributes_) {
            delete it.second;
        }
    }
}

void VariantContext::determineType() {
    if (this->type_ == VariantContextType::NULL_TYPE) {
        switch (getNAllele()) {
        case 0:
        {
            std::cerr << "[VariantContext::determineType] Error! request "
                << "type of VariantContext with no alleles!" << std::endl;
            std::exit(1);
        }
        case 1:
        {
            this->type_ = VariantContextType::NO_VARIATION;
            break;
        }
        default:
        {
            determinePolymorphicType();
            break;
        }
        }
    }
}


void VariantContext::determinePolymorphicType()
{
    type_ = VariantContextType::NULL_TYPE;
    Allele ref = getReference();

    const std::vector<Allele> &alleles = getAlleles();
    for (size_t i = 1; i < alleles.size(); ++i) {
        VariantContextType biallelicType = typeOfBiallelicVariant(
            ref, alleles[i]);
        if (type_ == VariantContextType::NULL_TYPE) {
            type_ = biallelicType;
        } else if (biallelicType != type_) {
            type_ = VariantContextType::MIXED;
            return;
        }
    }
}


VariantContextType VariantContext::typeOfBiallelicVariant(
    const Allele &ref, const Allele &alt)
{
    if (alt.isSymbolic()) {
        return VariantContextType::SYMBOLIC;
    } else if (ref.length() == alt.length()) {
        return alt.length() == 1 ?
            VariantContextType::SNP : VariantContextType::INDEL;
    } else {
        return VariantContextType::INDEL;
    }
}