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
#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <cerrno>
#include <cstddef>
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
#include "htslib/kstring.h"
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


static char *find_chrom_header_line(char *s)
{
    char *nl;
    if (strncmp(s, "#CHROM\t", 7) == 0) return s;
    else if ((nl = strstr(s, "\n#CHROM\t")) != NULL) return nl+1;
    else return NULL;
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
    bcf_hdr_t *merged_hrd1 = NULL;
    for (size_t i=0; i<headers.size(); ++i)
        merged_hrd1 = bcf_hdr_merge(merged_hrd1, headers[i]);
    
    int res = 0;
    
    // clear all samples
    res |= bcf_hdr_set_samples(merged_hrd1, NULL, 0) < 0;
    if (res) {
        std::cerr << "[bcf_hdr_merge_add_samples] Error! can not clear bcf "
            << "header samples" << std::endl;
        std::exit(1);
    }

    // not use keep_samples bit set; keep all samples
    if (merged_hrd1->keep_samples) free(merged_hrd1->keep_samples);
    merged_hrd1->keep_samples = NULL;
    
    kstring_t merged_hrd1_str = {0,0,0};
    kstring_t result_str = {0,0,0};

    if (bcf_hdr_format(merged_hrd1, 1, &merged_hrd1_str) < 0) {
        std::cerr << "[bcf_hdr_merge_add_samples] Error! can not convert "
            << "merged header 1 to string" << std::endl;
        ks_free(&merged_hrd1_str);
        ks_free(&result_str);
        bcf_hdr_destroy(merged_hrd1);
        std::exit(1);
    }

    bcf_hdr_t *result = bcf_hdr_init("w");
    bcf_hdr_set_version(result, bcf_hdr_get_version(merged_hrd1));

    if (merged_hrd1_str.l > 0) merged_hrd1_str.s[merged_hrd1_str.l-1] = '\t';
    res |= kputs("FORMAT\t", &merged_hrd1_str) < 0;
    char *p = find_chrom_header_line(merged_hrd1_str.s);
    int i = 0, end = 8;
    while ((p = strchr(p, '\t')) != 0 && i < end) ++i, ++p;
    if (i != end) {
        std::cerr << "[bcf_hdr_merge_add_samples] Error! Wrong number of "
            << "columns in header #CHROM line" << std::endl;
        ks_free(&merged_hrd1_str);
        ks_free(&result_str);
        bcf_hdr_destroy(merged_hrd1);
        std::exit(1);
    }

    res |= kputsn(merged_hrd1_str.s, p - merged_hrd1_str.s, &result_str) < 0;
    for (auto const &s: samples) {
        res |= kputc('\t', &result_str) < 0;
        res |= kputs(s.c_str(), &result_str) < 0;
    }

    while (result_str.l &&
        (!result_str.s[result_str.l-1] || result_str.s[result_str.l-1]=='\n') )
        result_str.l--; // kill trailing zeros and newlines
    res |= kputc('\n',&result_str) < 0;

    if (res) {
        std::cerr << "[bcf_hdr_merge_add_samples] Error! "
            << strerror(errno) << std::endl;
        ks_free(&merged_hrd1_str);
        ks_free(&result_str);
        bcf_hdr_destroy(merged_hrd1);
        std::exit(1);
    }

    if (bcf_hdr_parse(result, result_str.s) < 0) {
        bcf_hdr_destroy(result);
        result = NULL;
    }

    ks_free(&merged_hrd1_str);
    ks_free(&result_str);
    bcf_hdr_destroy(merged_hrd1);

    return result;
}
