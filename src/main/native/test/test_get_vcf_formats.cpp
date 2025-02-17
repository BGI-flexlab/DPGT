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
#include "htslib/vcf.h"
#include "vcf/vcf_reader.hpp"
#include <cstdint>
#include <iostream>
#include <vector>


void getPLs(bcf1_t *variant, bcf_hdr_t *hdr, std::vector<std::vector<int>> &PLs)
{
    std::cout << variant->rid << ":" << variant->pos + 1 << "\t";
    int ndst = 0; int32_t *dst = NULL;
    int n;
    int nsamples = bcf_hdr_nsamples(hdr);
    if ( (n = bcf_get_format_int32(hdr, variant, "PL", &dst, &ndst)) > 0 )
    {
        n /= nsamples;
        for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
            int32_t *ptr = dst + i*n;
            std::vector<int> PL;
            for (int j = 0; j < n; ++j) {
                if (ptr[j] == bcf_int32_vector_end) break;
                PL.push_back(ptr[j]);
            }
            if (!PL.empty()) {
                for (auto p: PL) std::cout << p << ",";
            } else {
                std::cout << "." << std::endl;
            }
            std::cout << "||";
            PLs.push_back(PL);
        }
    }
    std::cout << std::endl;
    free(dst);
}


int main(int argc, char **argv) {
    VcfReader reader(argv[1], false);

    bcf1_t *record = bcf_init1();
    bcf_hdr_t *header = reader.header();
    kstring_t line = {0, 0, NULL};
    while ( reader.Read(record) ) {
        vcf_format1(header, record, ks_clear(&line));
        printf("%s", line.s);
        std::vector<std::vector<int>> PLs;
        getPLs(record, header, PLs);
    }

    bcf_hdr_destroy(header);

    if (line.s) free(line.s);
} 

