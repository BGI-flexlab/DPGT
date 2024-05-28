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
#include "vcf/vcf_id_table.hpp"
#include "lix/lix.h"


enum class VcfReaderQueryStatus: uint8_t {
    NOT,
    SUCCESS,
    FAIL
};


/**
 * A htslib based vcf reader for reading vcf, bgzip vcf and bcf samenessly.
 */
class VcfReader {
public:
    VcfReader() = default;
    
    /**
     * Construct a vcf reader by openning a vcf/bcf file, 
     * optionally with an tbi/csi index file.
     * Can not be used with a lix index file.
     */
    VcfReader(const std::string &file_name, bool require_index=false) {
        Open(file_name, "", require_index, false);
    }

    /**
     * Construct a vcf reader by openning a vcf/bcf file,
     * optionally with an tbi/csi/lix index file.
     */
    VcfReader(const std::string &file_name, const std::string &index_name,
        bool require_index=false, bool use_lix=false)
    {
        Open(file_name, index_name, require_index, use_lix);
    }

    VcfReader(const VcfReader &reader) = delete;
    VcfReader(VcfReader &&reader) noexcept;

    VcfReader &operator==(const VcfReader &reader) = delete;
    VcfReader &operator==(VcfReader &&reader) noexcept;

    ~VcfReader() {
        Close();
    }

    void Open(const std::string &file_name, const std::string &index_name,
        bool require_index=false, bool use_lix=false);

    void Close() {
        if (header_) bcf_hdr_destroy(header_);
        if (bcf_idx_) hts_idx_destroy(bcf_idx_);
        if (tbx_idx_) tbx_destroy(tbx_idx_);
        if (lix_idx_) lix_destroy(lix_idx_);
        if (file_ptr_) hts_close(file_ptr_);
        if (itr_) hts_itr_destroy(itr_);
        if (lix_itr_) lix_iter_destroy(lix_itr_);
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

    void setIdTables(VcfKeyMaps *vcf_key_maps_) {
        id_tables_ = VcfIdTables(header(), vcf_key_maps_);
    }

    VcfIdTables *getIdTables() {
        return &id_tables_;
    }

private:
    htsFile *file_ptr_ = nullptr;
    tbx_t *tbx_idx_ = nullptr;
    hts_idx_t *bcf_idx_ = nullptr;
    lix_t *lix_idx_ = nullptr;
    bcf_hdr_t *header_ = nullptr;
    hts_itr_t *itr_ = nullptr;
    lix_iter_t *lix_itr_ = nullptr;
    bool streaming_ = true;  // streaming reader, not support random access
    kstring_t tmps_ = {0, 0, NULL};
    const std::vector<SimpleInterval> *intervals_ = nullptr; // target intervals
    std::vector<SimpleInterval>::const_iterator intervals_iter_;
    VcfIdTables id_tables_;

    VcfReaderQueryStatus query_status_ = VcfReaderQueryStatus::NOT;

    std::string getFileName() {
        if (file_ptr_ != nullptr) {
            return file_ptr_->fn;
        } else {
            return "";
        }
    }

    void OpenVcfWithHtsIdx(
        const std::string &file_name, const std::string &index_name,
        bool require_index=false);
    
    void OpenVcfWithLixIdx(
        const std::string &file_name, const std::string &index_name,
        bool require_index=false);
    
    VcfReader &QueryiWithHtsIdx(int32_t tid, int64_t start, int64_t end);
    VcfReader &QuerysWithHtsIdx(const std::string &chrom, int64_t start, int64_t end);
    VcfReader &QueryiWithLixIdx(int32_t tid, int64_t start, int64_t end);
    VcfReader &QuerysWithLixIdx(const std::string &chrom, int64_t start, int64_t end);

};

#endif  // VCF_READER_HPP