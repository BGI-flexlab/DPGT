#include "tools/combine_gvcfs.hpp"
#include "common/cached_ref_seq.hpp"
#include "common/interval.hpp"
#include "common/simple_interval.hpp"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "vcf/allele.hpp"
#include "vcf/multi_vcf_reader.hpp"
#include "vcf/variant_builder.hpp"
#include "vcf/variant_context.hpp"
#include "vcf/variant_context_utils.hpp"
#include "vcf/vcf_attribute.hpp"
#include "vcf/vcf_constants.hpp"
#include "vcf/vcf_reader.hpp"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <fstream>
#include "gperftools/profiler.h"



void CombineGVCFs::apply(std::vector<std::shared_ptr<VariantContext>> &variants,
    ReferenceContext &reference_context)
{
    // Check that the input variant contexts do not contain MNPs as these may
    // not be properly merged
    for (auto &v: variants) {
        if (VariantContextUtils::isUnmixedMnpIgnoringNonRef(v)) {
            std::cerr << "Combining gVCFs containing MNPs is not supported."
                << "find MNP at " << v->getContigName() << ":" << v->getStart()
                << std::endl;
            std::exit(1);
        }
    }

    if (!variants_overlap_current_merge_.empty()) {
        SimpleInterval last = !prev_pos_.IsNull() &&
            prev_pos_.tid ==
            variants_overlap_current_merge_.front()->getContig()
            ? prev_pos_ :
            VariantContextUtils::variantContextToInterval(
                variants_overlap_current_merge_.front());
        int end = last.tid == reference_context.getTid() ?
            reference_context.getStart() - 1 :
            maxVariantContextsEnd(variants_overlap_current_merge_.begin(),
                variants_overlap_current_merge_.end());
        createIntermediateVariants(SimpleInterval(last.tid, last.start, end));
    }

    // clear variants_overlap_current_merge_ if variants are out of
    // reference context
    for (auto itr = variants_overlap_current_merge_.begin();
        itr != variants_overlap_current_merge_.end(); )
    {
        if ((*itr)->getContig() != reference_context.getTid() ||
            (*itr)->getEnd() < reference_context.getStart())
        {
            itr = variants_overlap_current_merge_.erase(itr);
        } else {
            ++itr;
        }
    }

    mergeWithNewVCs(variants, reference_context);

    if (stored_reference_context_.isNull() ||
        stored_reference_context_.getTid() != reference_context.getTid() ||
        stored_reference_context_.getWindow().end <
        reference_context.getWindow().end)
    {
        stored_reference_context_ = reference_context;
    }
}

void CombineGVCFs::finalize() {
    MultiVariantGroupOnStart::finalize();
    if (!variants_overlap_current_merge_.empty()) {
        const SimpleInterval last_interval(
            variants_overlap_current_merge_.front()->getContig(),
            variants_overlap_current_merge_.front()->getStart(),
            maxVariantContextsEnd(variants_overlap_current_merge_.begin(),
                variants_overlap_current_merge_.end()));
        createIntermediateVariants(last_interval);
    }
}


void CombineGVCFs::createIntermediateVariants(
    const SimpleInterval &interval_to_close)
{
    resizeReferenceIfNeeded(interval_to_close);

    RelativeBitSet sites_to_stop = getIntermediateStopSites(
        interval_to_close, multiple_at_witch_to_break_bands_);
    
    for (auto &vc: variants_overlap_current_merge_) {
        // GATK gvcf record has two alleles(ref, <NON-REF>) for non-variant site,
        // if number of alleles greater than 2, then there must be at least one
        // alt allele
        if (vc->getNAllele() > 2) {
            for (int64_t i = interval_to_close.start;
                i <= vc->getEnd() && i <= interval_to_close.end; ++i)
            {
                sites_to_stop.set(i);
            }
        } else if (vc->getEnd() <= interval_to_close.end) {
            sites_to_stop.set(vc->getEnd());
        }
    }

    for (int64_t i = 0; i < sites_to_stop.size(); ++i) {
        if (sites_to_stop.testIndex(i)) {
            int64_t stop_loc = sites_to_stop.get(i);
            SimpleInterval loc(interval_to_close.tid, stop_loc, stop_loc);
            if (isWithinInterval(loc)) {
                std::string ref_bases = stored_reference_context_.getBases().
                    substr(
                    stop_loc-stored_reference_context_.getWindow().start, 2);
                endPreviousStates(loc, ref_bases, {}, true);
            }
        }
    }
}


RelativeBitSet CombineGVCFs::getIntermediateStopSites(
    const SimpleInterval &interval_to_close,
    int64_t break_bands_multiple)
{
    RelativeBitSet sites_to_stop(
        interval_to_close.start, interval_to_close.end);
    if (break_bands_multiple > 0) {
        int64_t block_end_pos = 
            ((interval_to_close.start + 1) / break_bands_multiple + 1) *
            break_bands_multiple - 1;
        for (;block_end_pos <= interval_to_close.end;
            block_end_pos += break_bands_multiple)
        {
            sites_to_stop.set(block_end_pos);
        }
    }
    return sites_to_stop;
}


void CombineGVCFs::mergeWithNewVCs(
    const std::vector<std::shared_ptr<VariantContext>> &variants,
    ReferenceContext &reference_context)
{
    if (!variants.empty()) {
        if (!okayToSkipThisSite(variants, reference_context)) {
            const SimpleInterval &loc = reference_context.getInterval();
            const int64_t close_pos = loc.start - 1;
            SimpleInterval close_loc(loc.tid, close_pos, close_pos);
            if (close_pos >= 0) {
                if (isWithinInterval(close_loc)) {
                    endPreviousStates(close_loc,
                        reference_context.getBases().substr(1),
                        variants, false);
                } else {
                    prev_pos_ = std::move(close_loc);
                    ref_after_prev_pos_ =
                        reference_context.getBases().substr(1)[0];
                }
            }
        }
    }
    for (auto const &v: variants) variants_overlap_current_merge_.push_back(v);
    boost::dynamic_bitset<> new_samples_indices = getSampleIndices(variants);
    sample_indices_ |= new_samples_indices;
}

bool CombineGVCFs::okayToSkipThisSite(
    const std::vector<std::shared_ptr<VariantContext>> &variants,
    const ReferenceContext &reference_context)
{
    boost::dynamic_bitset<> intersection = getSampleIndices(variants);
    intersection &= sample_indices_;
    return !prev_pos_.IsNull() &&
        reference_context.getStart() == prev_pos_.start + 1 &&
        intersection.count() == 0;
}

std::set<std::string> CombineGVCFs::getSamples(
    const std::vector<std::shared_ptr<VariantContext>> &variants)
{
    std::set<std::string> samples;
    for (auto const &v: variants) {
        char **samples1 = v->getSamples();
        for (int i = 0; i < v->getNSamples(); ++i) {
            samples.insert(samples1[i]);
        }
    }
    return samples;
}

boost::dynamic_bitset<> CombineGVCFs::getSampleIndices(
    const std::vector<std::shared_ptr<VariantContext>> &variants)
{
    boost::dynamic_bitset<> indices;
    if (variants.empty()) return indices;
    indices = boost::dynamic_bitset<>(variants.front()->getMergedNSamples());
    for (auto const &v: variants) {
        const std::vector<int> &indices1 = v->getSampleIndices();
        for (int i = 0; i < v->getNSamples(); ++i) {
            indices.set(indices1[i]);
        }
    }
    return indices;
}

static
bool containsAllSamples(const std::set<std::string> &samples, char **qsamples,
    int qsamples_size)
{
    for (int i = 0; i < qsamples_size; ++i) {
        if (samples.find(qsamples[i]) == samples.end()) return false;
    }
    return true;
}

static
void removeAllSampleIndices(boost::dynamic_bitset<> &samples_indices,
    const std::vector<int> &qsample_indices)
{
    for (auto &it: qsample_indices) {
        samples_indices.reset(it);
    }
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

static
void removeAllSamples(std::set<std::string> &samples, char **qsamples,
    int qsamples_size)
{
    for (int i = 0; i < qsamples_size; ++i) {
        auto itr = samples.find(qsamples[i]);
        if (itr != samples.end()) {
            samples.erase(itr);
        }
    }
}

void CombineGVCFs::endPreviousStates(const SimpleInterval &pos,
    const std::string &ref_bases,
    const std::vector<std::shared_ptr<VariantContext>> &variants,
    bool force_output_at_current_pos)
{
    boost::dynamic_bitset<> new_sample_indices = getSampleIndices(variants);

    char ref_base = ref_bases[0];
    char ref_next_base = force_output_at_current_pos ?
        (ref_bases.size() > 1 ? ref_bases[1] : 'N') : ref_base;
    
    std::vector<std::shared_ptr<VariantContext>> stopped_vcs;
    stopped_vcs.reserve(variants_overlap_current_merge_.size());

    for (auto itr = variants_overlap_current_merge_.begin();
        itr != variants_overlap_current_merge_.end(); )
    {
        std::shared_ptr<VariantContext> &vc = *itr;
        if (vc->getStart() <= pos.start || vc->getContig()!=pos.tid) {
            stopped_vcs.push_back(vc);

            if (vc->getEnd() == pos.start ||
                (variants.size() > 0 && !force_output_at_current_pos &&
                containsAllSampleIndices(
                    new_sample_indices, vc->getSampleIndices())))
            {
                removeAllSampleIndices(sample_indices_, vc->getSampleIndices());
                itr = variants_overlap_current_merge_.erase(itr);
            } else {
                ++itr;
            }
        } else {
            ++itr;
        }
    }

    if (!stopped_vcs.empty() && (prev_pos_.IsNull() || pos > prev_pos_))
    {
        SimpleInterval closing_spot = SimpleInterval(
            pos.tid, pos.start, pos.start);

        VariantStringWithPos merged_var;

        if (containsTrueAltAllele(stopped_vcs)) {
            merged_var = reference_confident_var_merger_.merge(
                stopped_vcs, closing_spot, ref_base, getMergedHeader(),
                getKeyMaps(), &tmp_var_ks_);
        } else {
            merged_var = referenceBlockMerge(stopped_vcs, closing_spot.start);
        }

        vcf_writer_->write(merged_var.var_str, merged_var.rid,
            merged_var.start, merged_var.end);
        
        prev_pos_ = std::move(closing_spot);
        ref_after_prev_pos_ = ref_next_base;
    }
}


bool CombineGVCFs::containsTrueAltAllele(
    const std::vector<std::shared_ptr<VariantContext>> &variants)
{
    for (auto const &vc: variants) {
        if (vc->getNAllele() > 2) return true;
    }
    return false;
}


VariantStringWithPos CombineGVCFs::referenceBlockMerge(
    std::vector<std::shared_ptr<VariantContext>> &variants, int64_t end)
{
    std::shared_ptr<VariantContext> &first = variants.front();

    Allele ref_allele;
    int64_t start;
    if (prev_pos_.IsNull() || prev_pos_.tid != first->getContig() ||
        first->getStart() >= prev_pos_.start + 1)
    {
        start = first->getStart();
        ref_allele = first->getReference();
    } else {
        start = prev_pos_.start + 1;
        ref_allele = Allele(std::string({ref_after_prev_pos_}), true);
    }
    std::vector<Allele> alleles_to_use = {ref_allele,
        Allele(Allele::NON_REF_ALLELE)};

    std::map<std::string, VcfAttributeBase *> shared_attributes;  // INFO
    
    if ( !use_bp_resolution_ && end != start ) {
        int32_t *end_val = (int32_t *)malloc(sizeof(int32_t));
        end_val[0] = end + 1;
        VcfAttribute<int32_t> *end_attr = new VcfAttribute<int32_t>(
            VCFConstants::END_KEY, BCF_HT_INT, 1, BCF_VL_FIXED, end_val);
        shared_attributes[VCFConstants::END_KEY] = end_attr;
    }

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
        first->getContig(), start, end, std::move(alleles_to_use));
    variant_builder.setAttributes(&shared_attributes).setGenotypes(
        &merged_genotypes);
    variant_builder.makeString(ks_clear(&tmp_var_ks_));

    for (auto &it: merged_genotypes) {
        delete it;
    }

    for (auto &it: shared_attributes) {
        delete it.second;
    }

    std::string var_str = ks_c_str(&tmp_var_ks_);
    return {var_str, variant_builder.tid(), variant_builder.start(),
        variant_builder.end()+1};
}


int main(int argc, char **argv) {
    ProfilerStart("t1.prof");

    if (argc != 5) {
        std::cerr << "usage: prog <vcfs.fofn> <ref.fa> <outfile> <region>" << std::endl;
        std::exit(1);
    }

    std::ifstream in(argv[1]);
    std::vector<std::string> vcfs;
    std::string line;
    while (std::getline(in, line)) {
        vcfs.push_back(line);
    }

    MultiVcfReader reader(vcfs, true);

    CachedRefSeq cached_ref_seq(argv[2]);

    int32_t tid = 0;
    hts_pos_t beg = 0;
    hts_pos_t end = 0;
    std::vector<SimpleInterval> intervals;
    const char *reg_ret = argv[4];
    while ((reg_ret = hts_parse_region(reg_ret, &tid, &beg, &end,
        (hts_name2id_f)bcf_hdr_name2id, reader.header(),
        HTS_PARSE_ONE_COORD|HTS_PARSE_LIST)) != NULL)
    {
        intervals.emplace_back(tid, beg, end-1);
    }

    reader.QueryIntervals(intervals);

    VcfWriter vcf_writer(argv[3]);

    CombineGVCFs combiner(&reader, &cached_ref_seq, &vcf_writer);

    combiner.run();

    ProfilerStop();
    return 0;
}
