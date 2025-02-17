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
#include "tools/vcf_header_combiner.hpp"
#include "common/check_jemalloc.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


int main(int argc, char **argv) {
    check_jemalloc();
    std::ifstream in(argv[1]);
    std::vector<std::string> vcfpaths;
    std::string line;
    while (std::getline(in, line)) {
        vcfpaths.push_back(line);
    }
    std::string outpath = "combine.header.vcf.gz";
    combineVCFHeaders(vcfpaths, outpath);
    std::cout << "end" << std::endl;

    std::cout << "end1" << std::endl;
}
