#include "tools/reference_confident_variant_merger.hpp"
#include "common/simple_matrix.hpp"
#include "genotyper/genotype_likelihood_calculator.hpp"
#include "htslib/vcf.h"
#include "vcf/allele.hpp"
#include "vcf/gatk_vcf_constants.hpp"
#include "vcf/variant_builder.hpp"
#include "vcf/variant_context.hpp"
#include "vcf/variant_context_utils.hpp"
#include "common/utils.hpp"
#include "vcf/vcf_attribute.hpp"
#include "vcf/vcf_constants.hpp"
#include "vcf/vcf_shared_attribute.hpp"
#include <algorithm>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>


static inline bool onlyContainsSpanDel(const std::vector<Allele> &alleles) {
    for (size_t i = 1; i < alleles.size(); ++i) {
        if (alleles[i].getDisplayString() != Allele::SPAN_DEL_STRING) {
            return false;
        }
    }
    return true;
}


VariantStringWithPos ReferenceConfidentVariantMerger::merge(
    std::vector<std::shared_ptr<VariantContext>> &vcs,
    const SimpleInterval &loc,
    char ref_base,
    bcf_hdr_t *merged_header,
    VcfKeyMaps *key_maps,
    bool remove_non_ref_allele,
    bool skip_span_del_only,
    int max_alt_alleles,
    kstring_t *out_var_ks)
{
    Allele ref_allele = determineRefAlleleGivenRefBase(vcs, loc, ref_base);

    // new allele map to old indices
    std::map<Allele, int *> new_allele_map_indices;

    for (size_t i = 0; i < vcs.size(); ++i) {
        std::shared_ptr<VariantContext> &vc = vcs[i];
        std::vector<Allele> new_alleles1;
        new_alleles1.reserve(vc->getNAllele() - 1);
        // if this variant not start at current loc then it must be a spaning
        // event(SPAN-DEL or ref block)
        const bool isSpanningEvent = loc.start != vc->getStart();
        if (isSpanningEvent) {
            new_alleles1 = relaceWithNoCallsAndDels(vc);
        } else {
            new_alleles1 = remapAlleles(vc, ref_allele);
        }
        // new_alleles1 are transformed from old alleles
        for (size_t j = 1; j < new_alleles1.size(); ++j) {
            auto itr = new_allele_map_indices.find(new_alleles1[j]);
            if (itr == new_allele_map_indices.end())
            {
                int *old_indices = (int *)malloc(vcs.size()*sizeof(int));
                // init old indices to -1
                for (size_t k = 0; k < vcs.size(); ++k) {
                    old_indices[k] = -1;
                }
                old_indices[i] = j;
                new_allele_map_indices[new_alleles1[j]] = old_indices;
            } else {
                if (itr->second[i] < 0) {
                    // the old index for new allele and i-th vc is not set
                    itr->second[i] = j;
                } else {
                    // the old index for new allele and i-th vc have been setted
                    // this new allele must be a duplicated allele for i-th vc
                    itr->second[i] = calculateBetterDupAlleleIndex(
                        vc, itr->second[i], j);
                }
            }
        }
    }

    // insert ref allele map to old indices, note that ref allele always map to
    // zeros
    int *old_ref_indices = (int *)calloc(vcs.size(), sizeof(int));
    new_allele_map_indices[ref_allele] = old_ref_indices;

    // remove NON_CALL allele from new allele map
    auto itr = new_allele_map_indices.find(Allele::NO_CALL);
    if (itr != new_allele_map_indices.end()) {
        free(itr->second);
        new_allele_map_indices.erase(itr);
    }

    if (remove_non_ref_allele) {
        itr = new_allele_map_indices.find(Allele::NON_REF_ALLELE);
        if (itr != new_allele_map_indices.end()) {
            free(itr->second);
            new_allele_map_indices.erase(itr);
        }
    }

    // non_ref_allele_indices for each variant context
    int *non_ref_allele_indices = (int *)malloc(vcs.size()*sizeof(int));
    for (size_t i = 0; i < vcs.size(); ++i) {
        non_ref_allele_indices[i] = vcs[i]->getNAllele() - 1;
    }

    // all new alleles
    std::vector<Allele> new_alleles;
    new_alleles.reserve(new_allele_map_indices.size());

    // get allele index map by reshape new_allele_map_indices
    // row: variant context index, column: new alleles index
    SimpleMatrix<int> allele_index_map(
        vcs.size(), new_allele_map_indices.size());
    int i = 0;
    for (auto &it: new_allele_map_indices) {
        new_alleles.push_back(it.first);
        for (size_t j = 0; j < vcs.size(); ++j) {
            allele_index_map(j, i) = it.second[j];
            if (allele_index_map(j, i) < 0) {
                // if new allele not match any old allele for this variant
                // map it to <NON-REF> allele index(last allele of vc)
                allele_index_map(j, i) = non_ref_allele_indices[j];
            }
        }
        ++i;
    }

    if (skip_span_del_only && new_alleles.size() == 2 && onlyContainsSpanDel(new_alleles)) {
        free(non_ref_allele_indices);
        for (auto &itr: new_allele_map_indices) {
            free(itr.second);
        }
        return {NULL, 0, 0, 0};
    } 

    // fliter alt alleles if it > max_alt_alleles
    if (max_alt_alleles > 0 && (int)(new_alleles.size() - 1) > max_alt_alleles)
    {
        new_alleles = calculateMostLikelyAlleles(new_alleles, allele_index_map, vcs, max_alt_alleles, !remove_non_ref_allele);

        // reconstruct allele_index_map
        allele_index_map = SimpleMatrix<int>(vcs.size(), new_alleles.size());
        for (size_t i = 0; i < new_alleles.size(); ++i) {
            int *old_allele_indices = new_allele_map_indices.at(new_alleles[i]);
            for (size_t j = 0; j < vcs.size(); ++j) {
                allele_index_map(j, i) = old_allele_indices[j];
                if (allele_index_map(j, i) < 0) {
                    // if new allele not match any old allele for this variant
                    // map it to <NON-REF> allele index(last allele of vc)
                    allele_index_map(j, i) = non_ref_allele_indices[j];
                }
            }
        }
    }

    free(non_ref_allele_indices);

    // free new_allele_map_indices
    for (auto &itr: new_allele_map_indices) {
        free(itr.second);
    }
    new_allele_map_indices.clear();
    
    // merged genotypes
    std::vector<FlatGenotype *> merged_genotypes(
        key_maps->sample_key_map.size(), nullptr);

    int32_t *depth = (int32_t *)calloc(1, sizeof(int32_t));
    std::unordered_set<std::string> rsIDs;

    for (size_t i = 0; i < vcs.size(); ++i) {
        std::vector<FlatGenotype *> merged_genotypes1 =
            mergeReferenceConfidenceGenotypes(vcs[i], allele_index_map.row(i),
                new_alleles.size());
        for (auto &g: merged_genotypes1) {
            merged_genotypes[g->getSampleIdx()] = g;
        }

        depth[0] += calculateVCDepth(vcs[i]);

        if (loc.start != vcs[i]->getStart()) continue;

        std::vector<std::string> vc_rsIDs;
        boost::split(vc_rsIDs, vcs[i]->getID(), boost::is_any_of(","));
        for (auto &id: vc_rsIDs) {
            if (id != ".") rsIDs.insert(id);
        }
    }

    // merge shared attributes(INFO)
    std::map<std::string, std::vector<VcfAttributeBase *>> vcs_attributes =
        getSharedAttributesOfVCs(vcs, loc);
    std::map<std::string, VcfAttributeBase *> merged_attributes =
        mergeSharedAttributes(vcs_attributes, allele_index_map);
    
    std::set<std::string> filters;
    bool saw_pass_sample = false;
    for (size_t i = 0; i < vcs.size(); ++i) {
        if (saw_pass_sample) continue;  // no get more filters if have seen pass 
        std::vector<std::string> vc_filters = vcs[i]->getFilter();
        if (vc_filters.empty() ||
            vc_filters[0] == VCFConstants::PASSES_FILTERS_v4) {
            saw_pass_sample = true;
        } else {
            for (auto &f: vc_filters) filters.insert(f);
        }
    }

    merged_attributes["DP"] = new VcfAttribute<int32_t>(
        "DP", BCF_HT_INT, 1, BCF_VL_FIXED, depth);

    std::string merged_rsID = rsIDs.empty() ? "." : boost::join(rsIDs, ",");

    std::string merged_filter;
    if (saw_pass_sample || filters.empty()) {
        merged_filter = ".";
    } else {
        merged_filter = boost::join(filters, ",");
    }

    // TODO build merged variant
    VariantBuilder variant_builder(merged_header, key_maps,
        loc.tid, loc.start, loc.start+new_alleles.front().length(), new_alleles);
    variant_builder.setID(std::move(merged_rsID))
        .setFilter(std::move(merged_filter))
        .setGenotypes(&merged_genotypes)
        .setAttributes(&merged_attributes);


    variant_builder.makeString(ks_clear(out_var_ks));

    for (auto &it: merged_genotypes) {
        if (it == nullptr) continue;
        // ugly, need remember what is newed and delete them here
        if (it->hasPL()) delete it->getPL();
        if (it->hasAD()) delete it->getAD();
        delete it;
    }

    for (auto &it: merged_attributes) {
        delete it.second;
    }

    return {out_var_ks, variant_builder.tid(), variant_builder.start(),
        variant_builder.end()};
}


Allele ReferenceConfidentVariantMerger::determineRefAlleleGivenRefBase(
    std::vector<std::shared_ptr<VariantContext>> &vcs,
    const SimpleInterval &loc, char ref_base)
{
    Allele ref_allele = VariantContextUtils::determineRefAllele(vcs, loc);

    if (ref_allele.isNull()) {
        // all input variants not starts at loc, for example,
        // allele vcs are SPAN-DELETION or NON-REF variant
        // create reference allele by ref_base
        return Allele(std::string({ref_base}), true);
    } else {
        return ref_allele;
    }
}


std::vector<Allele>
ReferenceConfidentVariantMerger::relaceWithNoCallsAndDels(
    std::shared_ptr<VariantContext> &vc)
{
    std::vector<Allele> result(vc->getNAllele());

    // no call the reference allele
    result[0] = Allele::NO_CALL;

    const std::vector<Allele> &alleles = vc->getAlleles();
    for (size_t i = 1; i < alleles.size(); ++i) {
        if (alleles[i].isNonRefAllele()) {
            result[i] = alleles[i];
        } else if (alleles[i].length() < vc->getReference().length()) {
            result[i] = Allele::SPAN_DEL;
        } else {
            result[i] = Allele::NO_CALL;
        }
    }
    return result;
}


std::vector<Allele> ReferenceConfidentVariantMerger::remapAlleles(
    std::shared_ptr<VariantContext> &vc, const Allele &ref_allele)
{
    Allele vc_ref = vc->getReference();
    const int extra_base_count = ref_allele.length() - vc_ref.length();
    if (extra_base_count < 0) {
        std::cerr << "[ReferenceConfidentVariantMerger::remapAlleles] Error! "
            << "the wrong reference was selected" << std::endl;
        std::exit(1);
    }

    std::string extra_bases =
        ref_allele.getBases().substr(vc_ref.length());

    std::vector<Allele> result;
    result.push_back(ref_allele);

    const std::vector<Allele> &alleles = vc->getAlleles();
    for (size_t i = 1; i < alleles.size(); ++i) {
        const Allele &a = alleles[i];
        if (a.isSymbolic()) {
            result.push_back(a);
        } else if (a.equals(Allele::SPAN_DEL)) {
            result.push_back(a);
        } else if (a.isCalled()) {
            if (extra_base_count == 0) {
                result.push_back(a);
            } else {
                std::string new_bases = a.getBases() + extra_bases;
                result.push_back(Allele(new_bases, false));
            }
        } else {
            result.push_back(a);
        }
    }

    return result;
}


std::vector<FlatGenotype *>
ReferenceConfidentVariantMerger::mergeReferenceConfidenceGenotypes(
    std::shared_ptr<VariantContext> &vc, int *allele_index_map,
    int num_new_alleles)
{
    const int max_ploidy_add1 = vc->getMaxPloidy()+1;

    // genotype index maps from new genotype index to old genotype index for
    // different ploidy
    int **genotype_index_maps_by_ploidy =
        (int **)malloc((max_ploidy_add1)*sizeof(int *));
    for (int i = 0; i < max_ploidy_add1; ++i)
        genotype_index_maps_by_ploidy[i] = nullptr;
    int *genotype_index_maps_by_ploidy_size =
        (int *)calloc(max_ploidy_add1, sizeof(int));
    
    // TODO: 什么情况下旧allele数目小于新allele数目？
    int max_allele_count = std::max(vc->getNAllele(), num_new_alleles);

    const std::vector<Allele> &alleles = vc->getAlleles();

    std::vector<FlatGenotype *> results;
    
    std::vector<FlatGenotype *> genotypes = vc->getFlatGenotypes();
    for (auto g: genotypes) {
        const int ploidy = g->getPloidy();

        FlatGenotype *new_genotype = new FlatGenotype(*g);

        if (g->hasPL()) {
            if (genotype_index_maps_by_ploidy[ploidy] == nullptr) {
                genotype_index_maps_by_ploidy[ploidy] =
                    createGenotypeLikelihoodCalculator(ploidy, max_allele_count,
                        calculators_).genotypeIndexMap(allele_index_map,
                        num_new_alleles, calculators_, genotype_index_maps_by_ploidy_size[ploidy]);
            }
            int *genotype_index_map = genotype_index_maps_by_ploidy[ploidy];
            
            VcfAttribute<int32_t> *new_pl = createPL(
                g->getPL(), genotype_index_map, genotype_index_maps_by_ploidy_size[ploidy]);
            new_genotype->setPL(new_pl);

            VcfAttribute<int32_t> *new_ad = nullptr;
            if (g->hasAD()) {
                new_ad = createAD(
                    g->getAD(), allele_index_map, num_new_alleles);
                new_genotype->setAD(new_ad);
            }
        } else if (excludeFromAnnotation(g, alleles)) {
            new_genotype->getGT()->setToNoCall();
        }

        // all set to no call, not support other method
        new_genotype->getGT()->setToNoCall();

        results.push_back(new_genotype);
    }

    free(genotype_index_maps_by_ploidy_size);
    // free genotype_index_maps_by_ploidy
    for (int i = 0; i < max_ploidy_add1; ++i) {
        if (genotype_index_maps_by_ploidy[i])
            free(genotype_index_maps_by_ploidy[i]);
    }
    free(genotype_index_maps_by_ploidy);

    return results;
}


VcfAttribute<int32_t> *ReferenceConfidentVariantMerger::createPL(
    VcfAttribute<int32_t> *old_pl, int *genotype_index_map, int size)
{
    return createNewAttributeFromOldAttribute(old_pl, genotype_index_map, size);
}


VcfAttribute<int32_t> *ReferenceConfidentVariantMerger::createAD(
    VcfAttribute<int32_t> *old_ad, int *allele_index_map, int size)
{
    return createNewAttributeFromOldAttribute(old_ad, allele_index_map, size);
}


bool ReferenceConfidentVariantMerger::excludeFromAnnotation(FlatGenotype *g,
    const std::vector<Allele> &var_alleles)
{
    return g->isHomRef(var_alleles) && g->hasPL() &&
        ((g->hasDP() && g->getDP()->operator[](0) == 0) || !g->hasDP()) &&
        g->hasGQ() && g->getGQ()->operator[](0) == 0;
}



int ReferenceConfidentVariantMerger::calculateBetterDupAlleleIndex(
    std::shared_ptr<VariantContext> &vc, int i, int j)
{
    const std::vector<FlatGenotype *> &genotypes = vc->getFlatGenotypes();
    int32_t min_pl_i = std::numeric_limits<int32_t>::max();
    int32_t min_pl_j = std::numeric_limits<int32_t>::max();
    for (auto g: genotypes) {
        if (!g->hasPL()) continue;
        auto calculator = createGenotypeLikelihoodCalculator(
            g->getPloidy(), vc->getNAllele(), calculators_);
        const int hom_i_gt_idx = calculator.allelesToIndex(
            {i, g->getPloidy()});
        const int hom_j_gt_idx = calculator.allelesToIndex(
            {j, g->getPloidy()});
        const int pl_i = (*g->getPL())[hom_i_gt_idx];
        const int pl_j = (*g->getPL())[hom_j_gt_idx];
        if (pl_i < min_pl_i) min_pl_i = pl_i;
        if (pl_j < min_pl_j) min_pl_j = pl_j;
    }
    return min_pl_i <= min_pl_j ? i : j;
}



int ReferenceConfidentVariantMerger::calculateVCDepth(
    std::shared_ptr<VariantContext> &vc)
{
    const std::map<std::string, VcfAttributeBase *> &shared_attributes =
        vc->getSharedAttributes();
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
    std::vector<FlatGenotype *> genotypes = vc->getFlatGenotypes();
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


std::map<std::string, std::vector<VcfAttributeBase *>>
ReferenceConfidentVariantMerger::getSharedAttributesOfVCs(
    std::vector<std::shared_ptr<VariantContext>> &vcs,
    const SimpleInterval &loc)
{
    std::map<std::string, std::vector<VcfAttributeBase *>> result;
    for (size_t i = 0; i < vcs.size(); ++i) {
        std::shared_ptr<VariantContext> &vc = vcs[i];
        // not add span-event attributes for merge
        if (vc->getStart() != loc.start) continue;
        const std::map<std::string, VcfAttributeBase *> &attributes =
            vc->getSharedAttributes();
        for (auto &a: attributes) {
            auto itr = result.find(a.first);
            if (itr == result.end()) {
                result[a.first] = std::vector<VcfAttributeBase *>(
                    vcs.size(), nullptr);
            }
            result[a.first][i] = a.second;
        }
    }
    return result;
}


std::map<std::string, VcfAttributeBase *>
ReferenceConfidentVariantMerger::mergeSharedAttributes(
    std::map<std::string, std::vector<VcfAttributeBase *>> &vcs_attributes,
    const SimpleMatrix<int> &allele_index_map)
{
    std::map<std::string, VcfAttributeBase *> result;
    for (auto &it: vcs_attributes) {
        if (VcfSharedAttributeConstants::STALE_KEYS.find(it.first) !=
            VcfSharedAttributeConstants::STALE_KEYS.end())
        {
            // skip 'stale' keys
            continue;
        }

        if (VcfSharedAttributeConstants::REDUCIBLE_KEYS.find(it.first) !=
            VcfSharedAttributeConstants::REDUCIBLE_KEYS.end())
        {
            // merge reducible attributes
            for (auto a: it.second) {
                if (a != nullptr) {
                   result[it.first] = a->merge(it.second, &allele_index_map);
                   break;
                }
            }
        } else {
            // merge other attributes by calculate median of attributes values
            VcfAttributeBase *first = nullptr;
            int16_t value_type;
            for (auto a: it.second) {
                if (a != nullptr)
                {
                    first = a;
                    value_type = a->type();
                    break;
                }
            }

            if (value_type == BCF_HT_INT) {
                std::vector<int32_t> values;
                for (auto a: it.second) {
                    if (a != nullptr) {
                        VcfAttribute<int32_t> *ac =
                            dynamic_cast<VcfAttribute<int32_t> *>(a);
                        for (int i = 0; i < ac->size(); ++i) {
                            values.push_back(ac->value()[i]);
                        }
                    }
                }
                int32_t *merged_value = (int32_t *)malloc(sizeof(int32_t));
                if (values.size() == 1) {
                    merged_value[0] = values[0];
                } else {
                    std::sort(values.begin(), values.end());
                    merged_value[0] = values[values.size()/2];  // get median
                }

                result[it.first] = new VcfAttribute<int32_t>(
                        first->key(), first->type(), 1, first->sizeType(),
                        merged_value);

            } else if (value_type == BCF_HT_REAL) {
                std::vector<float> values;
                for (auto a: it.second) {
                    if (a != nullptr) {
                        VcfAttribute<float> *ac =
                            dynamic_cast<VcfAttribute<float> *>(a);
                        for (int i = 0; i < ac->size(); ++i) {
                            values.push_back(ac->value()[i]);
                        }
                    }
                }
                float *merged_value = (float *)malloc(sizeof(float));
                if (values.size() == 1) {
                    merged_value[0] = values[0];
                } else {
                    std::sort(values.begin(), values.end());
                    merged_value[0] = values[values.size()/2];  // get median
                }

                result[it.first] = new VcfAttribute<float>(
                        first->key(), first->type(), 1, first->sizeType(),
                        merged_value);
            } // not merge if value_type = BCF_HT_STR
        }
    }
    return result;
}


std::vector<Allele> ReferenceConfidentVariantMerger::calculateMostLikelyAlleles(
    const std::vector<Allele> &new_alleles,
    const SimpleMatrix<int> &allele_index_map,
    const std::vector<std::shared_ptr<VariantContext>> &vcs,
    int num_alt_alleles_to_keep,
    bool has_non_ref_allele)
{
    int non_ref_allele_idx = -1;

    if (has_non_ref_allele) {
        for (int i = (int)new_alleles.size() - 1; i > -1; --i)
        {
            if (new_alleles[i].isNonRefAllele()) {
                non_ref_allele_idx = i;
                break;
            }
        }
    }

    if (non_ref_allele_idx > -1 &&
        static_cast<int>(new_alleles.size()) - 2 <= num_alt_alleles_to_keep)
    {
        return new_alleles;
    }

    std::vector<double> likelihood_sums = calculateLikelihoodSums(new_alleles, allele_index_map, vcs);
    auto cmp = [&likelihood_sums](int i, int j) -> bool { return likelihood_sums[i] > likelihood_sums[j]; };
    std::vector<int> indices(new_alleles.size(), 0);
    for (int i = 0; i < (int)new_alleles.size(); ++i) {
        indices[i] = i;
    }
    std::sort(indices.begin()+1, indices.end(), cmp);  // sork by likelihood sums, ref allele always at front
    std::sort(indices.begin()+1, indices.begin()+1+num_alt_alleles_to_keep);  // sort first num_alt_alleles_to_keep, to keep orders in new_alleles
    std::vector<Allele> results;
    results.reserve(num_alt_alleles_to_keep+2);
    int num_alt_alleles_to_keep_p1 = num_alt_alleles_to_keep + 1;
    bool has_add_non_ref = false;
    for (int i = 0; i < num_alt_alleles_to_keep_p1; ++i) {
        if (indices[i] == non_ref_allele_idx) {
            has_add_non_ref = true;
        }
        results.push_back(new_alleles[indices[i]]);
    }
    if (!has_add_non_ref && non_ref_allele_idx > -1) {
        results.push_back(new_alleles[non_ref_allele_idx]);
    }
    return results;
}


static int *reverseArrayIndexAndValue(const int *a, int l, int s, int *rev, int &rev_l) {
    for (int i = 0; i < rev_l; ++i) {
        rev[i] = -1;  // init values to -1
    }
    const int rev_l1 = rev_l - 1;
    int j;
    for (int i = 0; i < l; ++i) {
        j = a[i];
        if (j >= rev_l1) {
            const int new_size = j + 10;
            rev = (int *)realloc(rev, new_size*sizeof(int));
            for (int k = rev_l; k < new_size; ++k) {
                rev[k] = -1;  // init new allocated values to -1
            }
            rev_l = new_size;
        }
        if (j != s) rev[j] = i;  // keep rev[s] == -1
    }
    return rev;
}

std::vector<double> ReferenceConfidentVariantMerger::calculateLikelihoodSums(
    const std::vector<Allele> &new_alleles,
    const SimpleMatrix<int> &allele_index_map,
    const std::vector<std::shared_ptr<VariantContext>> &vcs)
{
    std::vector<double> likelihood_sums = std::vector<double>(new_alleles.size(), 0.0);
    int *old_to_new_index_map = (int *)malloc(new_alleles.size()*sizeof(int));
    int old_to_new_index_map_size = new_alleles.size();
    for (size_t i = 0; i < vcs.size(); ++i) {
        std::shared_ptr<VariantContext> vc = vcs[i];
        if (vc->getNAllele() < 3) {
            // assume vc is a gvcf variant record, number of allele < 3 means it is not a true variant
            continue;
        }
        old_to_new_index_map = reverseArrayIndexAndValue(allele_index_map.row(i), allele_index_map.cols(),
            vc->getNAllele()-1, old_to_new_index_map, old_to_new_index_map_size);
        for (auto &genotype: vc->getFlatGenotypes()) {
            GenotypeLikelihoods gls = genotype->getLikelihoods();
            if (gls.empty()) {
                continue;
            }
            const EigenArrayXd &log10_likelihoods = gls.log10_likelihoods();
            int index_of_max_gl;
            log10_likelihoods.maxCoeff(&index_of_max_gl);
            if (index_of_max_gl == 0) {
                continue;  // skip when hom-ref gt has max gl
            }
            double gl_diff_bt_ref_and_best = log10_likelihoods[index_of_max_gl] - log10_likelihoods[0];
            int ploidy = genotype->getPloidy() > 0 ? genotype->getPloidy() : default_ploidy_;

            // allele counts(genotype) of the max gl
            std::vector<int> allele_counts =
                createGenotypeLikelihoodCalculator(ploidy, vc->getNAllele(), calculators_)
                .genotypeAlleleCountsAt(index_of_max_gl)
                .alleleCountsByIndex(vc->getNAllele()-1);

            for (int a = 1; a < (int)allele_counts.size(); ++a) {
                if (allele_counts[a] > 0 && a < old_to_new_index_map_size && old_to_new_index_map[a] > 0) {
                    likelihood_sums[old_to_new_index_map[a]] += gl_diff_bt_ref_and_best;
                }
            }
        }
    }
    free(old_to_new_index_map);
    return likelihood_sums;
}

