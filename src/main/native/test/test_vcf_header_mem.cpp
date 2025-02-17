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
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

int help(const std::vector<std::string> &vcfpaths) {
    bcf_hdr_t **headers = (bcf_hdr_t **)malloc(vcfpaths.size()*sizeof(bcf_hdr_t *));
    char fmode[5];
    strcpy(fmode, "r");
    for (int i = 0; i < (int)vcfpaths.size(); ++i) {
        vcf_open_mode(fmode+1, vcfpaths[i].c_str(), NULL);
        htsFile *fp = hts_open(vcfpaths[i].c_str(), fmode);
        if ( ! fp ) {
            std::cerr << "[VCFHeaderCombiner_Combine] Error! Fail to open "
                << vcfpaths[i] << std::endl;
            std::exit(1);
        }
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        headers[i] = hdr;
        hts_close(fp);
    }

    for (int i = 0; i < (int)vcfpaths.size(); ++i) {
        bcf_hdr_destroy(headers[i]);
        // free(headers[i]);
        headers[i] = nullptr;
    }

    free(headers);

    std::cout << "held done" << std::endl;
}


int main(int args, char **argv) {
    std::ifstream in(argv[1]);
    std::vector<std::string> vcfpaths;
    std::string line;
    while (std::getline(in, line)) {
        vcfpaths.push_back(line);
    }
    help(vcfpaths);
    std::cout << "end" << std::endl;
}
