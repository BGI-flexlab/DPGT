#ifndef DPGT_GENOTYPE_LIKELIHOODS_HPP
#define DPGT_GENOTYPE_LIKELIHOODS_HPP

#include <limits>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

#include "boost/math/special_functions/binomial.hpp"   // binomial coefficient

#include "common/math_utils.hpp"
#include "Eigen/Core"
#include "Eigen/Dense"


class GenotypeLikelihoods;

/**
 * Cache of genotype number.
 */
class GenotypeNumLikelihoodsCache
{
private:
    static const int kDefaultNumberOfAlleles;
    static const int kDefaultPloidy;

    EigenArrayXXi static_cache;

    void Initialize(int num_alleles, int ploidy);

public:
    GenotypeNumLikelihoodsCache();
    
    /**
     * @param num_alleles max number of alleles for this cache. > 0
     * @param ploidy max ploidy for this cache. > 0
     */
    GenotypeNumLikelihoodsCache(int num_alleles, int ploidy);
    
    ~GenotypeNumLikelihoodsCache() {}

    /**
     * Get number of likelihoods given number of alleles and ploidy
     * @param num_alleles number of alleles. > 0 && <= max_alleles
     * @param ploidy  ploidy. > 0 && <= max_ploidy
     * @return the number of genotype likelihoods.
     */
    int GetNumLikelihoods(int num_alleles, int ploidy) const;

    int MaxCachedNumAlleles() const;

    int MaxCachedPloidy() const;
};


class PLIndexToAlleleIndicesCache {
public:
    using ArrayType =
        Eigen::Array<std::vector<std::vector<int>>, Eigen::Dynamic, Eigen::Dynamic>;

private:
    static const int kDefaultNumberOfAlleles;
    static const int kDefaultPloidy;

    ArrayType static_cache;

    void Initialize(int num_alleles, int ploidy);

public:
    PLIndexToAlleleIndicesCache();

    PLIndexToAlleleIndicesCache(int num_alleles, int ploidy);

    ~PLIndexToAlleleIndicesCache() {}

    std::vector<std::vector<int>> GetPLIndexToAlleleIndices(
        int num_alleles, int ploidy) const;
    
    int MaxCachedNumAlleles() const;

    int MaxCachedPloidy() const;
};


/**
 * Represent genotype likelihoods of arbitrary ploidy and allele number.
 */
class GenotypeLikelihoods
{
private:
    // The maximum number of diploid alternate alleles that we can represent as
    // genotype likelihoods.
    static const int kMaxDiploidAltAlleles;
    static GenotypeNumLikelihoodsCache num_likelihoods_cache;

    static PLIndexToAlleleIndicesCache pl_idx_to_allele_idxs_cache;

    EigenArrayXd log10_likelihoods_;

public:
    GenotypeLikelihoods() = default;
    GenotypeLikelihoods(const EigenArrayXd &log10_likelihoods);

    /**
     * Initialize from PL string. eg. 100,0,200
     * @param PL_string
     */
    GenotypeLikelihoods(const std::string &PL_string);

    static GenotypeLikelihoods fromPLs(const EigenArrayXi &PLs);

    static const int kMaxPL = std::numeric_limits<int>::max();
    static int CalcNumLikelihoods(int num_alleles, int ploidy);

    static int GetNumLikelihoods(int num_alleles, int ploidy);

    /**
     * Calculate PL index to alleles indices list(order of PLs) using
     * recursive method.
     * eg. for 1 alternative allele, diploid, the order of PLs is:
     * 00, 01, 11
     *     for 1 alternative allele, 3 ploid, the order of PLs is:
     * 000, 001, 011, 111
     * The implementation method is described in VCF4.3 specification document.
     * @param alt_alleles number of alternative alleles exclude the reference allele.
     * @param ploidy number of ploid.
     * @param pl_index_to_allele_indices result.
     */
    static void CalculatePLIndexToAlleleIndices(int alt_alleles, int ploidy,
        std::vector<std::vector<int>> &pl_index_to_allele_indices);
    
    static void GetPLIndexToAlleleIndices(int alt_alleles, int ploidy,
        std::vector<std::vector<int>> &pl_index_to_allele_indices);

    /**
     * Calculate PL index given 2 allele indices. Only works for diploid sample.
     * The implementation method is described in VCF4.3 specification document.
     * @param allele1_index index of allele1.
     * @param allele2_index index of allele2.
     * @return the PL index in the PL array.
     */
    static int CalculatePLIndex(int allele1_index, int allele2_index);

    /**
     * Calculate allele conditional probability given the genotype, thus,
     * Calculate P(a|G), a is allele, G is genotype.
     * @param allele_index allele index.
     * @param genotype genotype represent by allele indices. Not empty!
     * @return P(a|G).
     */
    static double CalculateAlleleProbOnGenotype(int allele_index,
            const std::vector<int> &genotype);

    /**
     * Return GLs given PLs.
     * @param PLs phred-scaled genotype likelihoods.
     * @return GLs.
     */
    static EigenArrayXd PLsToGLs(const EigenArrayXd &PLs);

    /**
     * Return adjusted PLs given GLs.
     * @param GLs log10 genotype likelihoods.
     * @return PLs.
     */
    static EigenArrayXi GLsToPLs(const EigenArrayXd &GLs);

    /**
     * Calculate adjusted PLs of double type given GLs.
     */
    static EigenArrayXd GLsToPLsDouble(const EigenArrayXd &GLs);

    static EigenArrayXi PLsDoubleToPLsInt(const EigenArrayXd &PLs);

    /**
     * @return the PL string.
     */
    std::string ToPLString() const;

    EigenArrayXd &log10_likelihoods();
    const EigenArrayXd &log10_likelihoods() const;

    bool empty() const {
        return log10_likelihoods_.size() == 0;
    }
};

#endif  // DPGT_GENOTYPE_LIKELIHOODS_HPP
