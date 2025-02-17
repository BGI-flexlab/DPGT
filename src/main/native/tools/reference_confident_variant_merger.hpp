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

    /**
     * Merge variants
     * @param vcs input variant contexts
     * @param loc locus
     * @param ref_base reference base start at locus
     * @param bcf_hdr_t merged header
     * @param key_maps vcf samples and format key maps for each vc,
     * key is sample name and format id, value is corresponding sample index
     * and format index in the header
     * @param remove_non_ref_allele if true then remove <NON_REF> allele from merged variant
     * @param skip_span_del_only if true and merged alleles only contains REF and SPAN_DEL(*),
     * then return a NULL variant string 
     * @param max_alt_alleles if max_alt_alleles is -1 then not filter alt alleles;
     * if max_alt_alleles > 0 then filter alt alleles to keep max_alt_alleles number of alt alleles
     * @param out_var_ks kstring_t buffer for write merged variant string, its data will be cleared before write
     * any new data
     * @return merged variant string and position(tid, start, end)
     */
    VariantStringWithPos merge(
        std::vector<std::shared_ptr<VariantContext>> &vcs,
        const SimpleInterval &loc,
        char ref_base,
        bcf_hdr_t *merged_header,
        VcfKeyMaps *key_maps,
        bool remove_non_ref_allele,
        bool skip_span_del_only,
        int max_alt_alleles,
        kstring_t *out_var_ks);

private:
    int default_ploidy_ = 2;
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

    /**
     * filter alleles by take top-N most likely alleles
     * @param new_alleles new alleles
     * @param allele_index_map new to old allele index map,
     * row: vc index, column: new allele index, value: old allele index
     * @param vcs variants to be merged
     * @param num_alt_alleles_to_keep number of alt alleles to be kept
     * @param has_non_ref_allele if <NON_REF> is in new_alleles
     */
    std::vector<Allele> calculateMostLikelyAlleles(
        const std::vector<Allele> &new_alleles,
        const SimpleMatrix<int> &allele_index_map,
        const std::vector<std::shared_ptr<VariantContext>> &vcs,
        int num_alt_alleles_to_keep,
        bool has_non_ref_allele);
    
    std::vector<double> calculateLikelihoodSums(
        const std::vector<Allele> &new_alleles,
        const SimpleMatrix<int> &allele_index_map,
        const std::vector<std::shared_ptr<VariantContext>> &vcs);
};



#endif  // DPGT_REFERENCE_CONFIDENT_VARIANT_MERGER
