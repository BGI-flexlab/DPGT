#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "htslib/vcf.h"
#include "vcf/hts_vcf_utils.hpp"



bcf_hdr_t *bcf_hdr_rename_sample(bcf_hdr_t *header,
    const std::vector<std::pair<int, std::string>> &old_idx_new_name_pairs)
{
    // not rename if input a empty list
    if (old_idx_new_name_pairs.empty()) return header;

    int nsamples = bcf_hdr_nsamples(header);
    std::vector<std::string> new_samples;
    for (int i = 0; i < nsamples; ++i) {
        new_samples.push_back(header->samples[i]);
    }

    for (auto const &p: old_idx_new_name_pairs) {
        if (p.first >= nsamples) {
            // out of range error
            std::cerr << "[bcf_hdr_rename_sample] Error! old index out of range."
                << "old header have " << nsamples << " samples, but input old"
                << " index " << p.first << std::endl;
            std::exit(1);
        }
        new_samples[p.first] = p.second;
    }

    // clear all samples
    int res = bcf_hdr_set_samples(header, NULL, 0);
    if (res < 0) {
        std::cerr << "[bcf_hdr_rename_sample] Error! can not clear bcf header "
            << "samples" << std::endl;
        std::exit(1);
    }

    // add new samples
    for (auto const &s: new_samples) {
        res = bcf_hdr_add_sample(header, s.c_str());
        if (res < 0) {
            std::cerr << "[bcf_hdr_rename_sample] Error! can not add sample "
                << "to bcf header " << std::endl;
            std::exit(1);
        }
    }

    return header;
}


void bcf_hdr_rename_dup_samples(std::vector<bcf_hdr_t *> headers,
    bool allow_dup_samples) {
    std::map<std::string, int> sample_counts;
    for (bcf_hdr_t *h: headers) {
        const int nsamples = bcf_hdr_nsamples(h);
        for (int i = 0; i < nsamples; ++i) {
            if (sample_counts.find(h->samples[i]) == sample_counts.end()) {
                sample_counts[h->samples[i]] = 1;
            } else {
                ++sample_counts[h->samples[i]];
            }
        }
    }

    std::map<std::string, int> dup_samples;
    for (auto const &it: sample_counts) {
        if (it.second > 1) dup_samples[it.first] = 0;
    }

    if (dup_samples.empty()) return;  // not find any duplicated samples

    // there are duplicated samples, but not allow duplicated samples
    // print error message and exit
    if (!allow_dup_samples) {
        std::cerr << "[bcf_hdr_rename_dup_samples] Error! find duplicated "
            "sample names in vcf files." << std::endl;
        std::exit(1);
    }

    for (bcf_hdr_t *h: headers) {
        const int nsamples = bcf_hdr_nsamples(h);
        std::vector<std::pair<int, std::string>> old_idx_new_name_pairs;
        for (int i = 0; i < nsamples; ++i) {
            if (dup_samples.find(h->samples[i]) != dup_samples.end()) {
                if (dup_samples[h->samples[i]] != 0) {
                    std::ostringstream ostring;
                    ostring << h->samples[i] << "-"
                        << dup_samples[h->samples[i]];
                    old_idx_new_name_pairs.push_back(
                        std::make_pair(i, ostring.str()));
                }
                ++dup_samples[h->samples[i]];
            }
        }
        bcf_hdr_rename_sample(h, old_idx_new_name_pairs);
    }
}


bcf_hdr_t *bcf_hdr_merge_add_samples(std::vector<bcf_hdr_t *> headers)
{
    std::set<std::string> samples;
    for (bcf_hdr_t *h: headers) {
        const int nsamples = bcf_hdr_nsamples(h);
        for (int i = 0; i < nsamples; ++i) {
            samples.insert(h->samples[i]);
        }
    }

    // merge header
    bcf_hdr_t *result = NULL;
    for (size_t i=0; i<headers.size(); ++i)
        result = bcf_hdr_merge(result, headers[i]);
    
    // clear all samples
    int res = bcf_hdr_set_samples(result, NULL, 0);
    if (res < 0) {
        std::cerr << "[bcf_hdr_merge_add_samples] Error! can not clear bcf "
            << "header samples" << std::endl;
        std::exit(1);
    }

    // not use keep_samples bit set; keep all samples
    if (result->keep_samples) free(result->keep_samples);
    result->keep_samples = NULL;

    // add new samples
    for (auto const &s: samples) {
        res = bcf_hdr_add_sample(result, s.c_str());
        if (res < 0) {
            std::cerr << "[bcf_hdr_merge_add_samples] Error! can not add sample"
                << " to bcf header " << std::endl;
            std::exit(1);
        }
        res = bcf_hdr_sync(result);
        if (res < 0) {
            std::cerr << "[bcf_hdr_merge_add_samples] Error! can not sync after"
                << "add sample to bcf header " << std::endl;
            std::exit(1);
        }
    }

    return result;
}


