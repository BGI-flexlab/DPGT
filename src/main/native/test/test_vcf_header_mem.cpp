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
