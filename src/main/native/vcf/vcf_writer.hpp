#ifndef DPGT_VCF_WRITER_HPP
#define DPGT_VCF_WRITER_HPP

#include <cstdint>
#include <iostream>
#include <string>
#include "htslib/vcf.h"


/**
 * @brief vcf/bcf wruter based on htslib
 */
class VcfWriter {
private:
    htsFile *fp_ = nullptr;
    std::string idx_fn_ = "";
    bcf_hdr_t *header_ = nullptr;

    bcf1_t *v_ = nullptr;
    kstring_t tmp_var_ks_ = {0, 0, NULL};

public:
    VcfWriter() = default;
    VcfWriter(const std::string &fn);
    VcfWriter(const VcfWriter &other) = delete;
    VcfWriter(VcfWriter &&other) noexcept = delete;
    VcfWriter &operator=(const VcfWriter &other) = delete;
    VcfWriter &operator=(VcfWriter &&other) noexcept = delete;
    ~VcfWriter() { close(); }

    void open(const std::string &fn);

    void close();

    void writeHeader(bcf_hdr_t *hdr);

    void write(kstring_t *v_str,
        int32_t rid, int64_t start, int64_t end);
    void write(bcf1_t *v);
};


#endif  // DPGT_VCF_WRITER_HPP
