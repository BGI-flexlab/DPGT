#ifndef DPGT_FLAT_GENOTYPE_HPP
#define DPGT_FLAT_GENOTYPE_HPP

#include "common/utils.hpp"
#include "genotyper/genotype.hpp"
#include "vcf/vcf_attribute.hpp"
#include "vcf/vcf_id_table.hpp"
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>


class FlatGenotype: public Genotype {
private:
    VcfAttributeGT *GT_ = nullptr;
    std::vector<VcfAttributeBase *> attributes_;
    VcfFormatKeyMap *fmt_key_map_;

public:
    FlatGenotype() = default;

    FlatGenotype(int sample_idx,
        VcfAttributeGT *GT,
        std::vector<VcfAttributeBase *> genotype_attributes,
        VcfFormatKeyMap *fmt_key_map): Genotype(sample_idx), GT_(GT),
        attributes_(std::move(genotype_attributes)), fmt_key_map_(fmt_key_map)
    {}

    FlatGenotype(const FlatGenotype &other): Genotype(other) {
        this->GT_ = other.GT_;
        this->attributes_ = other.attributes_;
        this->fmt_key_map_ = other.fmt_key_map_;
    }

    FlatGenotype(FlatGenotype &&other) noexcept : Genotype(std::move(other)) {
        this->GT_ = other.GT_;
        this->attributes_ = std::move(other.attributes_);
        this->fmt_key_map_ = other.fmt_key_map_;
    }

    FlatGenotype &operator=(const FlatGenotype &other) {
        if (this != &other) {
            Genotype::operator=(other);
            this->GT_ = other.GT_;
            this->attributes_ = other.attributes_;
            this->fmt_key_map_ = other.fmt_key_map_;
        }
        return *this;
    }

    FlatGenotype &operator=(FlatGenotype &&other) noexcept {
        if (this != &other) {
            Genotype::operator=(std::move(other));
            this->GT_ = other.GT_;
            this->attributes_ = std::move(other.attributes_);
            this->fmt_key_map_ = other.fmt_key_map_;
        }
        return *this;
    }

    ~FlatGenotype() override {}

    VcfAttributeGT *getGT() const override {
        return GT_;
    }

    VcfAttribute<int32_t> *getGQ() const override {
        if (fmt_key_map_->GQ_IDX < 0) return nullptr;
        return dynamic_cast<VcfAttribute<int32_t> *>(
            attributes_[fmt_key_map_->GQ_IDX]);
    }

    VcfAttribute<int32_t> *getDP() const override {
        if (fmt_key_map_->DP_IDX < 0) return nullptr;
        return dynamic_cast<VcfAttribute<int32_t> *>(
            attributes_[fmt_key_map_->DP_IDX]);
    }

    VcfAttribute<int32_t> *getAD() const override {
        if (fmt_key_map_->AD_IDX < 0) return nullptr;
        return dynamic_cast<VcfAttribute<int32_t> *>(
            attributes_[fmt_key_map_->AD_IDX]);
    }

    VcfAttribute<int32_t> *getPL() const override {
        if (fmt_key_map_->PL_IDX < 0) return nullptr;
        return dynamic_cast<VcfAttribute<int32_t> *>(
            attributes_[fmt_key_map_->PL_IDX]);
    }

    VcfAttribute<int32_t> *getMIN_DP() const override {
        if (fmt_key_map_->MIN_DP_IDX < 0) return nullptr;
        return dynamic_cast<VcfAttribute<int32_t> *>(
            attributes_[fmt_key_map_->MIN_DP_IDX]);
    }

    FlatGenotype &setGT(VcfAttributeGT *GT) {
        GT_ = GT;
        return *this;
    }

    FlatGenotype &setGQ(VcfAttribute<int32_t> *GQ) {
        if (fmt_key_map_->GQ_IDX < 0) {
            std::cerr << "[FlatGenotype] Error! Can not setGQ, GQ FORMAT ID"
                << " is not in merged vcf header!" << std::endl;
            std::exit(1);
        }
        attributes_[fmt_key_map_->GQ_IDX] = GQ;
        return *this;
    }

    FlatGenotype &setDP(VcfAttribute<int32_t> *DP) {
        if (fmt_key_map_->DP_IDX < 0) {
            std::cerr << "[FlatGenotype] Error! Can not setDP, DP FORMAT ID"
                << " is not in merged vcf header!" << std::endl;
            std::exit(1);
        }
        attributes_[fmt_key_map_->DP_IDX] = DP;
        return *this;
    }

    FlatGenotype &setAD(VcfAttribute<int32_t> *AD) {
        if (fmt_key_map_->AD_IDX < 0) {
            std::cerr << "[FlatGenotype] Error! Can not setAD, AD FORMAT ID"
                << " is not in merged vcf header!" << std::endl;
            std::exit(1);
        }
        attributes_[fmt_key_map_->AD_IDX] = AD;
        return *this;
    }

    FlatGenotype &setPL(VcfAttribute<int32_t> *PL) {
        if (fmt_key_map_->PL_IDX < 0) {
            std::cerr << "[FlatGenotype] Error! Can not setPL, PL FORMAT ID"
                << " is not in merged vcf header!" << std::endl;
            std::exit(1);
        }
        attributes_[fmt_key_map_->PL_IDX] = PL;
        return *this;
    }

    std::vector<int> getAvailableKeyIndices() const {
        std::vector<int> result;
        for (int i = 0 ; i < static_cast<int>(attributes_.size()); ++i) {
            if (attributes_[i]) result.push_back(i);
        }
        return result;
    }

    const std::vector<VcfAttributeBase *> &getAttributes() const {
        return attributes_;
    }

    void getString(const std::set<int> &format_key_indices,
        int ploidy, kstring_t *s) const;
};



#endif  // DPGT_FLAT_GENOTYPE_HPP
