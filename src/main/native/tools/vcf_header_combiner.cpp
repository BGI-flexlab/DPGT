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
#include <iostream>
#include "tools/vcf_header_combiner.hpp"
#include "htslib/vcf.h"
#include "vcf/hts_vcf_utils.hpp"
#include "vcf/vcf_writer.hpp"


void combineVCFHeaders(const std::vector<std::string> &vcfpaths, const std::string &outpath)
{
    const int chunk_size = 100;
    std::vector<bcf_hdr_t *> mheaders; // merged header for each chunk
    std::vector<bcf_hdr_t *> headers;
    char fmode[5];
    strcpy(fmode, "r");
    int n = 1;
    for (int i = 0; i < (int)vcfpaths.size(); ++i) {
        vcf_open_mode(fmode+1, vcfpaths[i].c_str(), NULL);
        htsFile *fp = hts_open(vcfpaths[i].c_str(), fmode);
        if ( ! fp ) {
            std::cerr << "[VCFHeaderCombiner_Combine] Error! Fail to open "
                << vcfpaths[i] << std::endl;
            std::exit(1);
        }
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        headers.push_back(hdr);
        hts_close(fp);
        if (n >= chunk_size) {
            bcf_hdr_t *mhdr = bcf_hdr_merge_add_samples(headers);
            for (auto it: headers) {
                bcf_hdr_destroy(it);
            }
            headers.clear();
            mheaders.push_back(mhdr);
            n = 1;
        }
        ++n;
    }

    if (headers.size() > 0) {
        bcf_hdr_t *mhdr = bcf_hdr_merge_add_samples(headers);
        mheaders.push_back(mhdr);
    }
    bcf_hdr_t *result = bcf_hdr_merge_add_samples(mheaders);
    VcfWriter writer{outpath};
    writer.writeHeader(result);

    for (auto it: headers) {
        bcf_hdr_destroy(it);
    }
    headers.clear();

    for (auto it: mheaders) {
        bcf_hdr_destroy(it);
    }
    mheaders.clear();

    bcf_hdr_destroy(result);
}