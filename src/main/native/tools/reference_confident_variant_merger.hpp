#ifndef DPGT_REFERENCE_CONFIDENT_VARIANT_MERGER
#define DPGT_REFERENCE_CONFIDENT_VARIANT_MERGER


#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include "common/simple_matrix.hpp"
#include "htslib/vcf.h"
#include "vcf/allele.hpp"
#include "vcf/variant_context.hpp"
#include "vcf/variant_builder.hpp"
#include "common/simple_interval.hpp"
#include "genotyper/genotype_likelihood_calculators.hpp"
#include "genotyper/genotype_likelihood_calculator.hpp"
#include "vcf/vcf_attribute.hpp"


class ReferenceConfidentVariantMerger {
public:

VariantStringWithPos merge(
    std::vector<std::shared_ptr<VariantContext>> &vcs,
    const SimpleInterval &loc, char ref_base, bcf_hdr_t *merged_header,
    VcfKeyMaps *key_maps, kstring_t *out_var_ks);

private:

GenotypeLikelihoodCalculators calculators_;

/**
 * @brief determine reference allele given reference base
 * 
 * @param vcs variants
 * @param loc position
 * @param ref_base reference base
 * @return Allele* new reference allele
 */
Allele determineRefAlleleGivenRefBase(
    std::vector<std::shared_ptr<VariantContext>> &vcs,
    const SimpleInterval &loc, char ref_base);

/**
 * @brief replace deletion allele with SPAN-DEL(*), <NON-REF> is kept as it is,
 * other alleles is replaced with NO-CALL.
 * Note that GATK replace all deletion alleles(length < ref_allele.length)
 * with SPAN-DEL, this may not be right if there are multiple deletion alleles
 * of different length. For example:
 * 1    12345   ATTATTA A,ATTA
 * at 1:12346 position, GATK will replace 'A' and 'ATTA' with SPAN-DEL,
 * but 'ATTA' is not a SPAN-DEL
 * 
 * @param vc 
 * @return std::vector<Allele> 
 */
std::vector<Allele> relaceWithNoCallsAndDels(
    std::shared_ptr<VariantContext> &vc);

std::vector<Allele> remapAlleles(
    std::shared_ptr<VariantContext> &vc, const Allele &ref_allele);


std::vector<FlatGenotype *>
mergeReferenceConfidenceGenotypes(
    std::shared_ptr<VariantContext> &vc, int *allele_index_map,
    int num_new_alleles);


template<typename T>
static
VcfAttribute<T> *createNewAttributeFromOldAttribute(VcfAttribute<T> *old,
    int *new_to_old_map, int size)
{
    T *new_values = (T *)malloc(size*sizeof(T));
    for (int i = 0; i < size; ++i) {
        new_values[i] = (*old)[new_to_old_map[i]];
    }

    VcfAttribute<T> *new_atrribute = new VcfAttribute<T>(
        old->key(), old->type(), size, old->sizeType(), new_values);
    
    return new_atrribute;
}


/**
 * @brief create new PL from old PL using genotype index map
 * 
 * @param old_pl old PL
 * @param genotype_index_map genotype index map, genotype index map is a array
 * of which offset is new PL index and value is old PL index
 * @param size size of genotype index map
 * @return VcfAttribute<int32_t> * new PL
 */
VcfAttribute<int32_t> *createPL(VcfAttribute<int32_t> *old_pl,
    int *genotype_index_map, int size);


/**
 * @brief create new AD from old AD using allele index map
 * 
 * @param old_ad old AD
 * @param allele_index_map allele index map, allele index map is a array of 
 * which offset is new AD index and value is old AD index
 * @param size size of allele index map
 * @return VcfAttribute<int32_t> * new AD
 */
VcfAttribute<int32_t> *createAD(VcfAttribute<int32_t> *old_ad,
    int *allele_index_map, int size);


bool excludeFromAnnotation(FlatGenotype *g,
    const std::vector<Allele> &var_alleles);


int calculateBetterDupAlleleIndex(std::shared_ptr<VariantContext> &vc,
    int i, int j);

/**
 * @brief calculate variant depth by getting from INFO directly or sum over sample
 * genotypes
 */
int calculateVCDepth(std::shared_ptr<VariantContext> &vc);

/**
 * @brief get shared attributes from variants, return a map of key to vector of
 * attributes for each input variant
 */
std::map<std::string, std::vector<VcfAttributeBase *>>
getSharedAttributesOfVCs(std::vector<std::shared_ptr<VariantContext>> &vcs,
    const SimpleInterval &loc);

std::map<std::string, VcfAttributeBase *>
mergeSharedAttributes(
    std::map<std::string, std::vector<VcfAttributeBase *>> &vcs_attributes,
    const SimpleMatrix<int> &allele_index_map);

};



#endif  // DPGT_REFERENCE_CONFIDENT_VARIANT_MERGER
