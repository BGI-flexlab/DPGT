#include <cstdlib>
#include <sstream>
#include <vector>

#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast/try_lexical_convert.hpp"
#include "common/math_utils.hpp"
#include "genotyper/genotype_likelihoods.hpp"


const int GenotypeNumLikelihoodsCache::kDefaultNumberOfAlleles = 5;
const int GenotypeNumLikelihoodsCache::kDefaultPloidy = 2;

void GenotypeNumLikelihoodsCache::Initialize(int num_alleles, int ploidy)
{
    static_cache = EigenArrayXXi(num_alleles, ploidy);
    for (int i = 0; i < num_alleles; ++i)
    {
        for (int j = 0; j < ploidy; ++j)
        {
            static_cache(i, j) = GenotypeLikelihoods::CalcNumLikelihoods(
                i+1, j+1);
        }
    }
}

GenotypeNumLikelihoodsCache::GenotypeNumLikelihoodsCache()
{
    Initialize(kDefaultNumberOfAlleles, kDefaultPloidy);
}

GenotypeNumLikelihoodsCache::GenotypeNumLikelihoodsCache(
    int num_alleles, int ploidy)
{
    Initialize(num_alleles, ploidy);
}


int GenotypeNumLikelihoodsCache::GetNumLikelihoods(
    int num_alleles, int ploidy) const
{
    return static_cache(num_alleles - 1, ploidy - 1);
}

int GenotypeNumLikelihoodsCache::MaxCachedNumAlleles() const
{
    return static_cache.rows();
}

int GenotypeNumLikelihoodsCache::MaxCachedPloidy() const
{
    return static_cache.cols();
}

const int PLIndexToAlleleIndicesCache::kDefaultNumberOfAlleles = 5;
const int PLIndexToAlleleIndicesCache::kDefaultPloidy = 2;

void PLIndexToAlleleIndicesCache::Initialize(int num_alleles, int ploidy)
{
    static_cache = ArrayType(num_alleles, ploidy);
    for (int i = 0; i < num_alleles; ++i)
    {
        for (int j = 0; j < ploidy; ++j)
        {
            GenotypeLikelihoods::CalculatePLIndexToAlleleIndices(i, j+1,
                static_cache(i, j));
        }
    }
}

PLIndexToAlleleIndicesCache::PLIndexToAlleleIndicesCache() {
    Initialize(kDefaultNumberOfAlleles, kDefaultPloidy);
}

PLIndexToAlleleIndicesCache::PLIndexToAlleleIndicesCache(
    int num_alleles, int ploidy)
{
    Initialize(num_alleles, ploidy);
}

std::vector<std::vector<int>>
PLIndexToAlleleIndicesCache::GetPLIndexToAlleleIndices(
        int num_alleles, int ploidy) const
{
    return static_cache(num_alleles - 1, ploidy - 1);
}

int PLIndexToAlleleIndicesCache::MaxCachedNumAlleles() const
{
    return static_cache.rows();
}

int PLIndexToAlleleIndicesCache::MaxCachedPloidy() const
{
    return static_cache.cols();
}

GenotypeNumLikelihoodsCache
GenotypeLikelihoods::num_likelihoods_cache = GenotypeNumLikelihoodsCache();

PLIndexToAlleleIndicesCache
GenotypeLikelihoods::pl_idx_to_allele_idxs_cache =
    PLIndexToAlleleIndicesCache();

const int GenotypeLikelihoods::kMaxDiploidAltAlleles = 10;

GenotypeLikelihoods::GenotypeLikelihoods(const EigenArrayXd &log10_likelihoods):
    log10_likelihoods_(log10_likelihoods) {}

GenotypeLikelihoods::GenotypeLikelihoods(const std::string &PL_string)
{
    std::vector<std::string> fields;
    boost::split(fields, PL_string, boost::is_any_of(","));
    log10_likelihoods_ = EigenArrayXd(fields.size());
    for (int i = 0; i < log10_likelihoods_.size(); ++i)
    {
        char *end_ptr;
        log10_likelihoods_[i] = strtol(fields[i].c_str(), &end_ptr, 10);
    }
    log10_likelihoods_ *= -0.1;
}

GenotypeLikelihoods GenotypeLikelihoods::fromPLs(const EigenArrayXi &PLs)
{
    EigenArrayXd log10_likelihoods(PLs.size());
    for (int i = 0; i < log10_likelihoods.size(); ++i) {
        log10_likelihoods[i] = -0.1 * PLs[i];
    }
    return GenotypeLikelihoods(log10_likelihoods);
}

int GenotypeLikelihoods::CalcNumLikelihoods(int num_alleles, int ploidy)
{
    return static_cast<int>(
        boost::math::binomial_coefficient<double>(
            ploidy + num_alleles - 1, ploidy));
}

int GenotypeLikelihoods::GetNumLikelihoods(int num_alleles, int ploidy)
{
    if (num_alleles <= num_likelihoods_cache.MaxCachedNumAlleles()
        && ploidy <= num_likelihoods_cache.MaxCachedPloidy())
    {
        return num_likelihoods_cache.GetNumLikelihoods(num_alleles, ploidy);
    } else {
        return CalcNumLikelihoods(num_alleles, ploidy);
    }
}

EigenArrayXd GenotypeLikelihoods::PLsToGLs(const EigenArrayXd &PLs)
{
    EigenArrayXd GLs = PLs / -10.0;
    return GLs;
}

EigenArrayXi GenotypeLikelihoods::GLsToPLs(const EigenArrayXd &GLs)
{
    EigenArrayXi PLs = EigenArrayXi(GLs.size());
    int max_idx;
    double adjust = GLs.maxCoeff(&max_idx);
    for (int i = 0; i < GLs.size(); ++i)
    {
        PLs[i] = static_cast<int>(
                std::round(std::min(-10.0*(GLs[i]-adjust),
                        static_cast<double>(kMaxPL))));
    }
    return PLs;
}


EigenArrayXd GenotypeLikelihoods::GLsToPLsDouble(const EigenArrayXd &GLs)
{
    EigenArrayXd PLs = EigenArrayXd(GLs.size());
    int max_idx;
    double adjust = GLs.maxCoeff(&max_idx);
    for (int i = 0; i < GLs.size(); ++i)
    {
        PLs[i] = std::min(-10.0*(GLs[i]-adjust), static_cast<double>(kMaxPL));
    }
    return PLs;
}


EigenArrayXi GenotypeLikelihoods::PLsDoubleToPLsInt(const EigenArrayXd &PLs)
{
    return VecRound<EigenArrayXd, EigenArrayXi>(PLs);
}


static void CalculatePLIndexToAlleleIndicesImplementation(int alt_alleles,
    int ploidy, std::vector<std::vector<int>> &pl_index_to_allele_indices,
    std::vector<int> genotype)
{
    for (int a = 0; a <= alt_alleles; ++a)
    {
        std::vector<int> gt{a};
        for (auto &it: genotype) gt.push_back(it);
        if (ploidy == 1)
        {
            pl_index_to_allele_indices.push_back(gt);
        } else if (ploidy > 1)
        {
            CalculatePLIndexToAlleleIndicesImplementation(
                    a, ploidy - 1, pl_index_to_allele_indices, gt);
        }
    }
}

void GenotypeLikelihoods::CalculatePLIndexToAlleleIndices(int alt_alleles,
    int ploidy, std::vector<std::vector<int>> &pl_index_to_allele_indices)
{
    CalculatePLIndexToAlleleIndicesImplementation(alt_alleles, ploidy,
            pl_index_to_allele_indices, std::vector<int>{});
}

void GenotypeLikelihoods::GetPLIndexToAlleleIndices(int alt_alleles, int ploidy,
    std::vector<std::vector<int>> &pl_index_to_allele_indices)
{
    int num_alleles = alt_alleles + 1;
    if (num_alleles <= pl_idx_to_allele_idxs_cache.MaxCachedNumAlleles()
        && ploidy <= pl_idx_to_allele_idxs_cache.MaxCachedPloidy())
    {
        pl_index_to_allele_indices =
            pl_idx_to_allele_idxs_cache.GetPLIndexToAlleleIndices(
                num_alleles, ploidy);
    } else {
        CalculatePLIndexToAlleleIndices(alt_alleles, ploidy,
            pl_index_to_allele_indices);
    }
}

int GenotypeLikelihoods::CalculatePLIndex(int allele1_index, int allele2_index)
{
    if (allele1_index > allele2_index)
    {
        std::swap(allele1_index, allele2_index);
    }
    return allele2_index * (allele2_index + 1) / 2 + allele1_index;
}


std::string GenotypeLikelihoods::ToPLString() const
{
    std::ostringstream osstream;
    EigenArrayXi PLs = GLsToPLs(log10_likelihoods_);
    for (int i = 0; i < PLs.size() - 1; ++i)
    {
        osstream << static_cast<int>(PLs[i]) << ",";
    }
    osstream << PLs[PLs.size() - 1];
    return osstream.str();
}

double GenotypeLikelihoods::CalculateAlleleProbOnGenotype(
        int allele_index, const std::vector<int> &genotype)
{
    int ploidy = genotype.size();
    int count = 0;
    for (auto const &a: genotype)
    {
        if (a == allele_index) ++count;
    }
    return static_cast<double>(count) / ploidy;
}


EigenArrayXd & GenotypeLikelihoods::log10_likelihoods()
{
    return log10_likelihoods_;
}

const EigenArrayXd &GenotypeLikelihoods::log10_likelihoods() const
{
    return log10_likelihoods_;
}


bool IsNullGenotype(const std::vector<int> &gt)
{
    if (gt[0] < 0 || gt[1] < 0) return true;
    return false;
}

bool IsHomRefGenotype(const std::vector<int> &gt)
{
    for (auto i: gt)
    {
        if (i != 0) return false;
    }
    return true;
}

bool IsHetGenotype(const std::vector<int> &gt)
{
    if (gt[0] == 0 && gt[1] > 0) return true;
    if (gt[0] > 0 && gt[1] == 0) return true;
    return false;
}

bool IsHomAltGenotype(const std::vector<int> &gt)
{
    if (gt[0] == gt[1] && gt[0] > 0 && gt[1] > 0) return true;
    return false;
}

bool IsHetAltGenotype(const std::vector<int> &gt)
{
    if (gt[0] != gt[1] && gt[0] > 0 && gt[1] > 0) return true;
    return false;
}

