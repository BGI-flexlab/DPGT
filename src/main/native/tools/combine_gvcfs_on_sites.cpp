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
#include "tools/combine_gvcfs_on_sites.hpp"
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cstddef>
#include <list>
#include <string>
#include <vector>


void CombineGVCFsOnSites::apply(std::shared_ptr<VariantContext> &&var)
{
    while (var->getStart() > cur_site_) {
        apply(cur_site_, variants_overlap_cur_site_);
        // goto next site
        cur_site_ = vs_set_->nextSite(cur_site_);
        // if next site is not exist, return
        if (cur_site_ < 0) return;
        // remove variants before next site
        for (auto itr = variants_overlap_cur_site_.begin();
            itr != variants_overlap_cur_site_.end(); )
        {
            if ((*itr)->getEnd() < cur_site_) {
                itr = variants_overlap_cur_site_.erase(itr);
            } else {
                ++itr;
            }
        }
    }

    if (var->getEnd() < cur_site_) return;

    variants_overlap_cur_site_.push_back(var);
}


static
bool containsAllSampleIndices(const boost::dynamic_bitset<> &samples_indices,
    const std::vector<int> &qsample_indices)
{
    for (auto &it: qsample_indices) {
        if (!samples_indices.test(it)) return false;
    }
    return true;
}


void CombineGVCFsOnSites::apply(
    int64_t site,
    std::list<std::shared_ptr<VariantContext>> &variants)
{
    // filter variants, if multiple variant of the same single vcf overlap with
    // this site, only keep the variant with largest start
    boost::dynamic_bitset<> samples_indices(getMergedSamplesNumber());
    std::list<std::shared_ptr<VariantContext>> variants_list;
    for (auto ritr = variants.rbegin(); ritr != variants.rend(); ++ritr) {
        if ( !containsAllSampleIndices(
            samples_indices, (*ritr)->getSampleIndices()) )
        {
            variants_list.push_front(*ritr);
        }
        for (auto i: (*ritr)->getSampleIndices()) samples_indices.set(i);
    }

    std::vector<std::shared_ptr<VariantContext>> variants_vec{
        variants_list.size(), nullptr};
    int i = 0;
    for (auto v: variants_list) {
        variants_vec[i++] = v;
    }

    VariantStringWithPos merged_var;
    if (containsTrueAltAllele(variants)) {
        std::string ref_seq = cached_ref_seq_->getSubsequenceAt(
            chrom_, site, site);
        merged_var = reference_confident_var_merger_.merge(
            variants_vec, SimpleInterval{tid_, site, site}, ref_seq.front(),
            getMergedHeader(), getKeyMaps(), false, false, -1, &tmp_var_ks_);
    } else {
        merged_var = referenceBlockMerge(variants_vec, site);
    }

    if (merged_var.var_ks && merged_var.var_ks->l > 1) {
        vcf_writer_->write(merged_var.var_ks, merged_var.rid,
            merged_var.start, merged_var.end);
    }
}


void CombineGVCFsOnSites::finalize() {
    while (cur_site_ > 0) {
        apply(cur_site_, variants_overlap_cur_site_);
        cur_site_ = vs_set_->nextSite(cur_site_);
        // remove variants before next site
        for (auto itr = variants_overlap_cur_site_.begin();
            itr != variants_overlap_cur_site_.end(); )
        {
            if ((*itr)->getEnd() < cur_site_) {
                itr = variants_overlap_cur_site_.erase(itr);
            } else {
                ++itr;
            }
        }
    }
}


void CombineGVCFsOnSites::run() {
    if (cur_site_ < 0) return;  // return if input variant site bitset is empty
    std::shared_ptr<VariantContext> variant_context;
    while ((variant_context = reader_->Read()) != nullptr) {
        apply(std::move(variant_context));
        if (cur_site_ < 0) break;  // if no variant site left, break
    }
    finalize();
}


VariantStringWithPos CombineGVCFsOnSites::referenceBlockMerge(
    std::vector<std::shared_ptr<VariantContext>> &variants, int64_t site)
{
    std::string ref_seq = cached_ref_seq_->getSubsequenceAt(chrom_, site, site);
    Allele ref_allele = Allele(ref_seq, true);
    std::vector<Allele> alleles_to_use = {ref_allele,
        Allele(Allele::NON_REF_ALLELE)};
    
    std::map<std::string, VcfAttributeBase *> shared_attributes;  // INFO
    int32_t *end_val = (int32_t *)malloc(sizeof(int32_t));
    end_val[0] = site + 1;
    VcfAttribute<int32_t> *end_attr = new VcfAttribute<int32_t>(
            VCFConstants::END_KEY, BCF_HT_INT, 1, BCF_VL_FIXED, end_val);
    shared_attributes[VCFConstants::END_KEY] = end_attr;

    std::vector<FlatGenotype *> merged_genotypes(
        getMergedSamplesNumber(), nullptr);

    for (auto &v: variants) {
        for (auto &g: v->getFlatGenotypes()) {
            FlatGenotype *new_g = new FlatGenotype(*g);
            new_g->getGT()->setToNoCall();
            merged_genotypes[new_g->getSampleIdx()] = new_g;
        }
    }

    VariantBuilder variant_builder(getMergedHeader(), getKeyMaps(),
        tid_, site, site+1, std::move(alleles_to_use));
    variant_builder.setAttributes(&shared_attributes).setGenotypes(
        &merged_genotypes);
    variant_builder.makeString(ks_clear(&tmp_var_ks_));

    for (auto &it: merged_genotypes) {
        delete it;
    }

    for (auto &it: shared_attributes) {
        delete it.second;
    }

    return {&tmp_var_ks_, variant_builder.tid(), variant_builder.start(),
        variant_builder.end()};
}
