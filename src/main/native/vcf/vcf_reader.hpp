#ifndef VCF_READER_HPP
#define VCF_READER_HPP

#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <string>
#include <cinttypes>
#include <vector>

#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "common/simple_interval.hpp"


/**
 * A htslib based vcf reader for reading vcf, bgzip vcf and bcf samenessly.
 */
class VcfReader {
public:
    VcfReader() = default;
    VcfReader(const std::string &file_name, bool require_index=false) {
        Open(file_name, require_index);
    }

    VcfReader(const VcfReader &reader) = delete;
    VcfReader(VcfReader &&reader) noexcept;

    VcfReader &operator==(const VcfReader &reader) = delete;
    VcfReader &operator==(VcfReader &&reader) noexcept;

    ~VcfReader() {
        Close();
    }

    void Open(const std::string &file_name, bool require_index=false);

    void Close() {
        if (header_) bcf_hdr_destroy(header_);
        if (bcf_idx_) hts_idx_destroy(bcf_idx_);
        if (tbx_idx_) tbx_destroy(tbx_idx_);
        if (file_ptr_) hts_close(file_ptr_);
        if (itr_) hts_itr_destroy(itr_);
        if (tmps_.s) free(tmps_.s);
    }


    VcfReader &Queryi(int32_t tid, int64_t start, int64_t end);

    VcfReader &Querys(const std::string &chrom, int64_t start, int64_t end);

    VcfReader &QueryIntervals(const std::vector<SimpleInterval> &intervals);

    /**
     * Read one vcf record
     */
    bcf1_t *Read(bcf1_t *record);

    bcf_hdr_t *header() {
        return header_;
    }

    const std::vector<SimpleInterval> *getIntervals() const {
        return intervals_;
    }

private:
    htsFile *file_ptr_ = nullptr;
    tbx_t *tbx_idx_ = nullptr;
    hts_idx_t *bcf_idx_ = nullptr;
    bcf_hdr_t *header_ = nullptr;
    hts_itr_t *itr_ = nullptr;
    bool streaming_ = true;  // streaming reader, not support random access
    kstring_t tmps_ = {0, 0, NULL};
    const std::vector<SimpleInterval> *intervals_ = nullptr; // target intervals
    std::vector<SimpleInterval>::const_iterator intervals_iter_;
};

#endif  // VCF_READER_HPP