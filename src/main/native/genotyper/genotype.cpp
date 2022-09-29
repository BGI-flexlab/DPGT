#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include "genotype.hpp"
#include "vcf/allele.hpp"
#include "vcf/vcf_attribute.hpp"
#include "vcf/vcf_constants.hpp"


const std::vector<std::string> Genotype::PRIMARY_KEYS = {
    VCFConstants::GENOTYPE_FILTER_KEY,
    VCFConstants::GENOTYPE_KEY,
    VCFConstants::GENOTYPE_QUALITY_KEY,
    VCFConstants::DEPTH_KEY,
    VCFConstants::GENOTYPE_ALLELE_DEPTHS,
    VCFConstants::GENOTYPE_PL_KEY
};

const std::string Genotype::PHASED_ALLELE_SEPARATOR = "|";
const std::string Genotype::UNPHASED_ALLELE_SEPARATOR = "/";

GenotypeType Genotype::determineType(
    const std::vector<Allele> &var_alleles) const 
{
    VcfAttributeGT *gt = getGT();
    if (gt == nullptr) return GenotypeType::UNAVAILABLE;

    bool saw_no_call = false;
    bool saw_multiple_alleles = false;
    Allele first_call_allele;
    int first_call_allele_index = -1;

    int allele_index = 0;
    for (int i = 0; i < gt->size(); ++i) {
        allele_index = (*gt)[i];
        if (allele_index < 0) {
            saw_no_call = true;
            continue;
        }
        const Allele &allele = var_alleles[(*gt)[i]];
        if (allele.isNoCall()) {
            saw_no_call = true;
        } else if (first_call_allele.isNull()) {
            first_call_allele = allele;
            first_call_allele_index = (*gt)[i];
        } else if ((*gt)[i] != first_call_allele_index) {
            saw_multiple_alleles = true;
        }
    }

    if (saw_no_call) {
        if (first_call_allele.isNull()) return GenotypeType::NO_CALL;
        return GenotypeType::MIXED;
    }

    if (first_call_allele.isNull()) {
        std::cerr << "BUG: there are no alleles present in this genotype"
            << " but the alleles list is not null" << std::endl;
        std::exit(1);
    }

    if (saw_multiple_alleles) {
        return GenotypeType::HET;
    } else {
        if (first_call_allele.isReference()) {
            return GenotypeType::HOM_REF;
        } else {
            return GenotypeType::HOM_VAR;
        }
    }        
}
