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