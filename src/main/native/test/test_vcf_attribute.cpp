#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <set>
#include "htslib/vcf.h"
#include "vcf/vcf_attribute.hpp"
#include "vcf/vcf_reader.hpp"


std::vector<std::map<std::string, VcfAttributeBase *>>
createGenotypeAttributes(bcf1_t *variant_, bcf_hdr_t *header_)
{
    const int nsamples = bcf_hdr_nsamples(header_);
    std::vector<std::map<std::string, VcfAttributeBase *>> result(nsamples);

    // loop through all genotype fields
    for (int i = 0; i < variant_->n_fmt; ++i) {
        const char *key =
            header_->id[BCF_DT_ID][variant_->d.fmt[i].id].key;
        const int val_type = variant_->d.fmt[i].type;
        switch (val_type) {
            case BCF_BT_INT8:
            case BCF_BT_INT16:
            case BCF_BT_INT32:
            case BCF_BT_INT64:
            {
                // dst , ndst
                int ndst = 0;
                int32_t *dst = nullptr;
                int n = bcf_get_format_int32(
                    header_, variant_, key, &dst, &ndst);
                if (n > 0) {
                    n /= nsamples;
                    for (int j = 0; j < nsamples; ++j) {
                        int32_t *ptr = dst + j*n;
                        int k;
                        for (k = 0; k < n; ++k) {
                            if (ptr[k] == bcf_int32_vector_end) break;
                        }
                        const size_t s = k*sizeof(int32_t);
                        int32_t *tmp = (int32_t *)malloc(s);
                        memcpy(tmp, ptr, s);
                        const uint8_t l = bcf_hdr_id2length(
                            header_, BCF_HL_FMT, variant_->d.fmt[i].id);
                        if (strcmp(key, "GT") == 0) {
                            // GT coded phase into int32_t, so decoded it 
                            // here
                            result[j][key] = new VcfAttributeGT(
                                key, BCF_HT_INT, k, l, tmp);
                        } else {
                            result[j][key] = new VcfAttribute<int32_t>(
                                key, BCF_HT_INT, k, l, tmp);
                        }
                    }
                }
                free(dst);
                break;
            }
            case BCF_BT_FLOAT:
            {
                // dst , ndst
                int ndst = 0;
                float *dst = nullptr;
                int n = bcf_get_format_float(
                    header_, variant_, key, &dst, &ndst);
                if (n > 0) {
                    n /= nsamples;
                    for (int j = 0; j < nsamples; ++j) {
                        float *ptr = dst + j*n;
                        int k;
                        for (k = 0; k < n; ++k) {
                            if (ptr[k] == bcf_float_vector_end) break;
                        }
                        const size_t s = k*sizeof(float);
                        float *tmp = (float *)malloc(s);
                        memcpy(tmp, ptr, s);
                        const uint8_t l = bcf_hdr_id2length(
                            header_, BCF_HL_FMT, variant_->d.fmt[i].id);
                        result[j][key] = new VcfAttribute<float>(
                            key, BCF_HT_REAL, k, l, tmp);
                    }
                }
                free(dst);
                break;
            }
            case BCF_BT_CHAR:
            {
                // dst , ndst
                int ndst = 0;
                char **dst = nullptr;
                int n = bcf_get_format_string(
                    header_, variant_, key, &dst, &ndst);
                if (n > 0) {
                    for (int j = 0; j < nsamples; ++j) {
                        const uint8_t l = bcf_hdr_id2length(
                            header_, BCF_HL_FMT, variant_->d.fmt[i].id);
                        result[j][key] = new VcfAttributeString(
                            key, BCF_HT_STR, 1, l, dst[j]);
                    }
                }
                free(dst[0]); free(dst);
                break;
            }
            default:
            {
                std::cerr << "Error! Invalid format value type for key="
                    << key << " value type=" << val_type << std::endl;
                std::exit(1);
            }
        }
    }
    return result;
}


int main(int argc, char **argv) {
    // std::vector<VcfAttributeBase *> attributes;

    // attributes.push_back(new VcfAttribute<int>("FLAG", 0, BCF_HT_FLAG, nullptr));
    
    // int *dp = (int *)malloc(sizeof(int));
    // dp[0] = 36;
    // attributes.push_back(new VcfAttribute<int>("DP", 1, BCF_HT_INT, dp));

    // float *fv = (float *)malloc(sizeof(float)*2);
    // fv[0] = 1.0;
    // fv[1] = 2.0;
    // attributes.push_back(new VcfAttribute<float>("FV", 2, BCF_HT_REAL, fv));

    // const char **aux = (const char **)malloc(sizeof(const char *)*3);
    // aux[0] = "A";
    // aux[1] = "B";
    // aux[2] = "C";
    // attributes.push_back(new VcfAttribute<const char *>("AUX", 3, BCF_HT_STR, aux));

    // for (auto a: attributes) {
    //     std::cout << a->key() << "=" << a->getValueStr() << std::endl;
    //     delete a;
    // }

    VcfReader reader(argv[1], false);

    bcf1_t *record = bcf_init1();
    bcf_hdr_t *header = reader.header();
    kstring_t line = {0, 0, NULL};
    while ( reader.Read(record) ) {
        vcf_format1(header, record, ks_clear(&line));
        printf("%s", line.s);
        auto fmts = createGenotypeAttributes(record, header);
        for (auto f: fmts) {
            for (auto kv: f) {
                std::cout << kv.first << "=" << kv.second->getValueStr() << ";";
                delete kv.second;
            }
            std::cout << "\t";
        }
        std::cout << std::endl;
    }

    if (line.s) free(line.s);
    bcf_destroy1(record);

}

