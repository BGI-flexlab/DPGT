#ifndef DPGT_GENOTYPE_HPP
#define DPGT_GENOTYPE_HPP

#include <cstddef>
#include <unordered_map>
#include <vector>
#include <string>
#include <cstdint>
#include <set>


#include "common/math_utils.hpp"
#include "vcf/allele.hpp"
#include "genotyper/genotype_likelihoods.hpp"
#include "vcf/vcf_attribute.hpp"
#include "vcf/vcf_constants.hpp"



enum class GenotypeType:uint8_t {
    // null genotype type
    NULL_TYPE,
    // the sample is no-called(all alleles are NO_CALL)
    NO_CALL,
    // the sample is homozygous reference
    HOM_REF,
    // the sample is heterozygous, with at least one ref and at least of one alt in any order
    HET,
    // all alleles are non-reference
    HOM_VAR,
    // there is no allele data for this sample(alleles is empty)
    UNAVAILABLE,
    // same chromosomes are NO_CALL and others are called(eg. ./1)
    MIXED
};


/**
 * @brief This class encompasses all the basic information about a genotype.
 * It is immutable.
 */
class Genotype {
public:
    static const std::vector<std::string> PRIMARY_KEYS;

    static const std::string PHASED_ALLELE_SEPARATOR;
    static const std::string UNPHASED_ALLELE_SEPARATOR;

private:
    int sample_idx_ = -1;
    GenotypeType type_ = GenotypeType::NULL_TYPE;

public:
    Genotype() = default;

    Genotype(int sample_idx): sample_idx_(sample_idx) {}
    
    Genotype(const Genotype &other) {
        this->sample_idx_ = other.sample_idx_;
        this->type_ = other.type_;
    }

    Genotype(Genotype &&other) noexcept {
        this->sample_idx_ = other.sample_idx_;
        this->type_ = other.type_;
    }

    Genotype &operator=(Genotype &&other) noexcept {
        if (&other != this) {
            this->sample_idx_ = other.sample_idx_;
            this->type_ = other.type_;
        }
        return *this;
    }

    Genotype &operator=(const Genotype &other) {
        if (&other != this) {
            this->sample_idx_ = other.sample_idx_;
            this->type_ = other.type_;
        }
        return *this;
    }
    
    virtual ~Genotype() {}

    bool operator==(const Genotype &other) const {
        return compareTo(other) == 0;
    }

    bool operator>(const Genotype &other) const {
        return compareTo(other) > 0;
    }

    bool operator<(const Genotype &other) const {
        return compareTo(other) < 0;
    }

    bool operator>=(const Genotype &other) const {
        return compareTo(other) >= 0;
    }

    bool operator<=(const Genotype &other) const {
        return compareTo(other) <= 0;
    }

    int getSampleIdx() const {
        return sample_idx_;
    }
    
    /**
     * @brief Get the Alleles
     * 
     * @return std::vector<Allele> vector of allele for this genotype
     */
    std::vector<Allele> getAlleles(
        const std::vector<Allele> &var_alleles) const
    {
        std::vector<Allele> gt_alleles;
        VcfAttributeGT *gt = getGT();
        for (int i = 0; i < gt->size(); ++i) {
            gt_alleles.push_back(var_alleles[(*gt)[i]]);
        }
        return gt_alleles;
    }

    /**
     * @brief count allele times in this genotype by its index
     * 
     * @param allele_index input allele index
     * @return int number of times of input allele in this genotype
     */
    int countAlleleByIndex(int allele_index) {
        VcfAttributeGT *gt = getGT();
        int c = 0;
        for (int i = 0; i < gt->size(); ++i) {
            if ((*gt)[i] == allele_index) ++c;
        }
        return c;
    }

    /**
     * @brief Get the Allele by index, note that index is offset in ploidy, not
     * allele index.
     * For example:
     * for genotype: 1/2, given index 0 will get allele for allele index 1
     * given index 1 will get allele for allele index 2
     * 
     * @param i ploidy offset
     * @return Allele
     */
    const Allele &getAllele(int i,
        const std::vector<Allele> &var_alleles) const {
        VcfAttributeGT *gt = getGT();
        return var_alleles[(*gt)[i]];
    }
    
    /**
     * @return true if this genotype is phased(separated by "|")
     */
    bool isPhased() const {
        return getGT()->isPhased();
    }

    /**
     * @return int ploidy of this genotype
     */
    int getPloidy() const {
        return getGT()->size();
    }

    virtual VcfAttributeGT *getGT() const = 0;
    virtual VcfAttribute<int32_t> *getGQ() const = 0;
    virtual VcfAttribute<int32_t> *getDP() const = 0;
    virtual VcfAttribute<int32_t> *getAD() const = 0;
    virtual VcfAttribute<int32_t> *getPL() const = 0;
    virtual VcfAttribute<int32_t> *getMIN_DP() const = 0;

    bool hasGT() const {
        return getGT() != nullptr;
    }

    bool hasDP() const {
        return getDP() != nullptr;
    }

    bool hasAD() const {
        return getAD() != nullptr;
    }

    bool hasGQ() const {
        return getGQ() != nullptr;
    }

    bool hasPL() const {
        return getPL() != nullptr;
    }

    bool hasMIN_DP() const {
        return getMIN_DP() != nullptr;
    }

    /**
     * @brief Get genotype type, type is cached
     * @return GenotypeType 
     */
    GenotypeType getType(const std::vector<Allele> &var_alleles) {
        if (type_ != GenotypeType::NULL_TYPE) {
            return type_;
        }
        type_ = determineType(var_alleles);
        return type_;
    }

    /**
     * @return true if all alleles of this genotype are the same
     */
    bool isHom(const std::vector<Allele> &var_alleles) {
        return isHomRef(var_alleles) || isHomVar(var_alleles);
    }

    /**
     * @return true if all alleles of this genotype are the reference
     */
    bool isHomRef(const std::vector<Allele> &var_alleles) {
        return getType(var_alleles) == GenotypeType::HOM_REF;
    }

    /**
     * @return true if all alleles of this genotype are the same and
     * are not reference
     */
    bool isHomVar(const std::vector<Allele> &var_alleles) {
        return getType(var_alleles) == GenotypeType::HOM_VAR;
    }

    /**
     * @return true if there is allele that differs from other alleles of this
     * genotype
     */
    bool isHet(const std::vector<Allele> &var_alleles) {
        return getType(var_alleles) == GenotypeType::HET;
    }

    /**
     * @return true if all allele are not reference and there is allele that
     * differs from other alleles
     * if ploidy is less than 2 return false
     */
    bool isHetNonRef(const std::vector<Allele> &var_alleles) {
        if (getPloidy() < 2) return false;
        return getType(var_alleles) == GenotypeType::HET &&
            getAllele(0, var_alleles).isNonReference() &&
            getAllele(1, var_alleles).isNonReference();
    }

    bool isNoCall(const std::vector<Allele> &var_alleles) {
        return getType(var_alleles) == GenotypeType::NO_CALL;
    }

    bool isMisex(const std::vector<Allele> &var_alleles) {
        return getType(var_alleles) == GenotypeType::MIXED;
    }

    bool isAvailable(const std::vector<Allele> &var_alleles) {
        return getType(var_alleles) != GenotypeType::UNAVAILABLE;
    }

    // genotype likelihood related routines
    
    /**
     * @return true if this genotype has PL field
     */
    bool hasLikelihoods() const {
        return hasPL();
    }

    std::string getLikelihoodString() const {
        std::ostringstream osstream;
        VcfAttribute<int32_t> *PL = getPL();
        const int last = PL->size() - 1;
        for (int i = 0; i < last; ++i)
        {
            osstream << static_cast<int>((*PL)[i]) << ",";
        }
        osstream << (*PL)[last];
        return osstream.str();
    }

    GenotypeLikelihoods getLikelihoods() const {
        VcfAttribute<int32_t> *PL = getPL();
        // convert to eigen array
        EigenArrayXi pl_array(PL->size());
        for (int i = 0; i < PL->size(); ++i) {
            pl_array[i] = (*PL)[i];
        }
        return hasLikelihoods() ? GenotypeLikelihoods::fromPLs(pl_array) :
            GenotypeLikelihoods();
    }

    /**
     * @brief are PLs non-informative? if PLs are empty or all fields of PLs
     * are 0, then it is non-informative
     * @return true if PLs are empty or all fields of PLs are 0.
     */
    bool isNonInformative() const {
        if (!hasPL()) {
            return true;
        } else {
            VcfAttribute<int32_t> *PL = getPL();
            for (int i = 0; i < PL->size(); ++i) {
                if ((*PL)[i] != 0) return false;
            }
            return true;
        }
    }

    int compareTo(const Genotype &other) const {
        return Compare(this->sample_idx_, other.sample_idx_);
    }

protected:
    GenotypeType determineType(
        const std::vector<Allele> &var_alleles) const;


};



#endif // DPGT_GENOTYPE_HPP
