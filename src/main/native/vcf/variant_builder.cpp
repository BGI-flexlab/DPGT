#include <cstddef>
#include <ios>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include <set>
#include <sstream>
#include <iomanip>
#include "boost/algorithm/string/join.hpp"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "vcf/hts_vcf_utils.hpp"
#include "vcf/variant_builder.hpp"
#include "vcf/vcf_constants.hpp"
#include "common/utils.hpp"



const float VariantBuilder::FLOAT_MISSING = std::numeric_limits<float>::min();

int VariantBuilder::getMaxPloidy(int default_ploidy)
{
    int max_ploidy = -1;
    std::vector<FlatGenotype *> &genotypes = getGenotypes();
    for (auto &it: genotypes)
    {
        if (it == nullptr) continue;
        if (it->getPloidy() > max_ploidy) max_ploidy = it->getPloidy();
    }
    if (max_ploidy == -1) return default_ploidy;
    return max_ploidy;
}

VariantBuilder::VariantBuilder(
    bcf_hdr_t *header, VcfKeyMaps *key_maps,
    int32_t tid, int64_t start, int64_t end, std::vector<Allele> alleles):
    header_(header), key_maps_(key_maps),
    tid_(tid), start_(start), end_(end), alleles_(std::move(alleles))
{
    qual_ = FLOAT_MISSING;
    filter_ = VCFConstants::MISSING_VALUE_v4;
}


VariantBuilder &VariantBuilder::setAttributes(
    std::map<std::string, VcfAttributeBase *> *attributes)
{
    attributes_ = attributes;
    return *this;
}

VariantBuilder &VariantBuilder::setGenotypes(
    std::vector<FlatGenotype *> *genotypes)
{
    genotypes_ = genotypes;
    return *this;
}


void VariantBuilder::makeString(kstring_t *out_ks) {
    kputs(bcf_hdr_id2name(header_, tid_), out_ks);
    kputc('\t', out_ks);
    kputl(start_ + 1, out_ks);
    kputc('\t', out_ks);
    kputs(id_.c_str(), out_ks);
    kputc('\t', out_ks);
    kputs(alleles_.front().getDisplayString().c_str(), out_ks);
    kputc('\t', out_ks);
    
    int last = alleles_.size() - 1;
    for (int i = 1; i < last; ++i) {
        kputs(alleles_[i].getDisplayString().c_str(), out_ks);
        kputc(',', out_ks);
    }
    kputs(alleles_[last].getDisplayString().c_str(), out_ks);
    kputc('\t', out_ks);

    kputs(filter_.c_str(), out_ks);
    kputc('\t', out_ks);

    if (qual_ == FLOAT_MISSING)
    {
        kputs(".\t", out_ks);
    } else {
        kputd(qual_, out_ks);
        kputc('\t', out_ks);
    }

    std::map<std::string, VcfAttributeBase *> &shared_attributes =
        getAttributes();

    if (!shared_attributes.empty()) {
        auto last_itr = shared_attributes.end();
        --last_itr;
        for (auto i = shared_attributes.begin(); i != last_itr; ++i) {
            kputs(i->first.c_str(), out_ks);
            kputc('=', out_ks);
            i->second->getValueStr(out_ks);
            kputc(';', out_ks);
        }
        kputs(last_itr->first.c_str(), out_ks);
        kputc('=', out_ks);
        last_itr->second->getValueStr(out_ks);
        kputc('\t', out_ks);
    } else {
        kputs(".\t", out_ks);
    }

    std::vector<FlatGenotype *> &genotypes = getGenotypes();

    std::set<int> format_key_indices;  // all format keys other than 'GT'
    for (auto &it: genotypes) {
        if (it == nullptr) continue;
        for (int i = 0; i < (int)it->getAttributes().size(); ++i)
        {
            if (it->getAttributes()[i]) format_key_indices.insert(i);
        }
    }

    // vcf specification required that 'GT' is the first field of FORMAT and 
    // FORMAT must have 'GT' in it
    kputs("GT", out_ks);
    if (!format_key_indices.empty()) {
        for (auto &i: format_key_indices) {
            kputc(':', out_ks);
            kputs(key_maps_->format_key_map.keys[i].c_str(), out_ks);
        }
    }

    kputc('\t', out_ks);

    int ploidy = getMaxPloidy();
    std::string sample;
    last = bcf_hdr_nsamples(header_) - 1;
    for (int i = 0; i < last; ++i) {
        if (genotypes[i] != nullptr) {
            genotypes[i]->getString(format_key_indices, ploidy, out_ks);
        } else {
            Utils::krepeatPuts(out_ks, VCFConstants::MISSING_VALUE_v4,
                VCFConstants::UNPHASED, ploidy);
        }
        kputc('\t', out_ks);
    }

    if (genotypes[last] != nullptr) {
        genotypes[last]->getString(format_key_indices, ploidy, out_ks);
    } else {
        Utils::krepeatPuts(out_ks, VCFConstants::MISSING_VALUE_v4,
            VCFConstants::UNPHASED, ploidy);
    }

    kputc('\n', out_ks);
}
