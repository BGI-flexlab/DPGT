#ifndef MULTI_VCF_READER_HPP
#define MULTI_VCF_READER_HPP

#include <cstdint>
#include <stack>
#include <vector>
#include <string>
#include <queue>
#include <memory>
#include "htslib/vcf.h"
#include "vcf/vcf_reader.hpp"
#include "vcf/vcf_ibuffer.hpp"
#include "vcf/variant_context.hpp"
#include "vcf/vcf_id_table.hpp"


/**
 * Read one variant per time, in position order, from multiple vcf files
 */
class MultiVcfReader {
public:
    MultiVcfReader() = default;
    MultiVcfReader(std::vector<std::string> files, bool require_index=false);

    MultiVcfReader(const MultiVcfReader &multi_vcf_reader) = delete;
    MultiVcfReader(MultiVcfReader &&multi_vcf_reader) noexcept;

    MultiVcfReader &operator==(const MultiVcfReader &multi_vcf_reader) = delete;
    MultiVcfReader &operator==(MultiVcfReader &&multi_vcf_reader) noexcept;

    ~MultiVcfReader();

    std::shared_ptr<VariantContext> First();

    std::shared_ptr<VariantContext> Read();

    bool IsEOF();

    void Queryi(int32_t tid, int64_t start, int64_t end);

    void Querys(const std::string &chrom, int64_t start, int64_t end);

    void QueryIntervals(const std::vector<SimpleInterval> &intervals);

    bcf_hdr_t *header();

    VcfKeyMaps *KeyMaps() {
        return vcf_key_maps_;
    }

    const std::vector<SimpleInterval> *getIntervals() const {
        if (readers_.empty()) return nullptr;
        return readers_.front()->getIntervals();
    }

private:
    std::vector<std::string> files_;
    std::vector<VcfReader *> readers_;
    std::priority_queue<VcfIBuffer *,
        std::vector<VcfIBuffer *>, VcfIBufferGreater> buffers_;
    std::vector<VcfIBuffer *> buffers_vec_;  // use for manage memory of buffers
    bcf_hdr_t *header_ = nullptr;
    VcfKeyMaps *vcf_key_maps_ = nullptr;
    
    void ReloadBuffers();

    void Free() {
        for (VcfReader *r: readers_) delete r;
        for (VcfIBuffer *b: buffers_vec_) delete b;
        bcf_hdr_destroy(header_);
        if (vcf_key_maps_) delete vcf_key_maps_;
    }
};



#endif  // MULTI_VCF_READER_HPP
