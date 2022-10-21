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
