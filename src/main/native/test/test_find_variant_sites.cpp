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
#include <fstream>
#include <string>
#include <vector>
#include "variant/variant_site_set.hpp"
#include "variant_site_finder/variant_site_finder.hpp"


int main(int argc, char **argv) {
    std::ifstream input = std::ifstream(argv[1]);
    std::vector<std::string> vcf_files;
    std::string line;
    while(std::getline(input, line)) {
        vcf_files.push_back(line);
    }

    VariantSiteSet vset = FindVariantSite(vcf_files, "chr21", 0, 46709982);

    return 0;
}