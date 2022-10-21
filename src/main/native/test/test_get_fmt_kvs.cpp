#include <map>
#include <set>
#include <string>
#include <vector>
#include "htslib/vcf.h"
#include "vcf/vcf_reader.hpp"



std::vector<std::map<std::string, std::string>> createGenotypeAttributes(
    bcf1_t *variant_, bcf_hdr_t *header_,
    const std::set<std::string> &exclude_keys)
{
    const int nsamples = bcf_hdr_nsamples(header_);
    std::vector<std::map<std::string, std::string>> result(nsamples);
    int n = 0;
    char **values = nullptr;
    // loop through all genotype fields
    for (int i = 0; i < variant_->n_fmt; ++i) {
        const char *key =
            header_->id[BCF_DT_ID][variant_->d.fmt[i].id].key;
        // key not in exclude keys set
        if (exclude_keys.find(key) == exclude_keys.end()) {
           if ( bcf_get_format_string(header_, variant_, key, &values, &n) > 0)
           {
               for (int j = 0; j < nsamples; ++j) {
                   result[j][key] = values[i];
               }
           }
        }
    }
    return result;
}


int main(int argc, char **argv) {
    VcfReader reader(argv[1], false);

    bcf1_t *record = bcf_init1();
    bcf_hdr_t *header = reader.header();
    kstring_t line = {0, 0, NULL};
    while ( reader.Read(record) ) {
        vcf_format1(header, record, ks_clear(&line));
        printf("%s", line.s);
        auto fmts = createGenotypeAttributes(record, header, {});
        for (auto f: fmts) {
            for (auto kv: f) {
                std::cout << kv.first << "=" << kv.second << ";";
            }
            std::cout << "\t";
        }
        std::cout << std::endl;
    }

    bcf_hdr_destroy(header);

    if (line.s) free(line.s);
}

