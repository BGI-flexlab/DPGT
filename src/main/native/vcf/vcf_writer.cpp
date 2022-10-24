#include "vcf/vcf_writer.hpp"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include <boost/algorithm/string/predicate.hpp>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>


VcfWriter::VcfWriter(const std::string &fn)
{
    open(fn);
}


void VcfWriter::open(const std::string &fn) {
    if (boost::ends_with(fn, ".vcf.gz")) {
        fp_ = bcf_open(fn.c_str(), "w1z");
        idx_fn_ = fn + ".tbi";
        // hts_set_threads(fp_, 1);
    } else if (boost::ends_with(fn, ".vcf")) {
        fp_ = bcf_open(fn.c_str(), "w");
    } else if (boost::ends_with(fn, ".bcf")) {
        fp_ = bcf_open(fn.c_str(), "wb1");
        idx_fn_ = fn + ".csi";
        // hts_set_threads(fp_, 1);
    }
}


void VcfWriter::close()
{
    if (!idx_fn_.empty())
    {
        if (bcf_idx_save(fp_) < 0) {
            fprintf(stderr, "[VcfWriter] Error saving index\n");
            std::exit(1);
        }
    }
    if (fp_) hts_close(fp_);
    if (v_) bcf_destroy1(v_);
    ks_free(&tmp_var_ks_);
}

void VcfWriter::writeHeader(bcf_hdr_t *hdr)
{
    if (bcf_hdr_write(fp_, hdr) < 0) {
        std::cerr << "[VcfWriter::writeHeader] Failed to write vcf header."
            << std::endl;
        std::exit(1);
    }
    header_ = hdr;
    if (!idx_fn_.empty()) {
        // init vcf index for write index on the fly
        if (bcf_idx_init(fp_, hdr, 0, idx_fn_.c_str()) < 0) {
            std::cerr << "[VcfWriter::writeHeader] Failed to init vcf index."
                << std::endl;
            std::exit(1);
        }
    }
}


void VcfWriter::write(const std::string &v_str, int32_t rid, int64_t start,
    int64_t end)
{
    if (fp_->format.format == vcf || fp_->format.format == text_format) {
        int ret;
        if (fp_->format.compression != no_compression) {
            ret = bgzf_write(fp_->fp.bgzf, v_str.c_str(), v_str.size());
        }
        else {
            ret = hwrite(fp_->fp.hfile, v_str.c_str(), v_str.size());
        }
        
        if (ret != static_cast<int>(v_str.size())) {
            std::cerr << "[VcfWriter::write] IO Error!" << std::endl;
            std::exit(1);
        }

        if (fp_->idx) {
            int tid;
            if ((tid = hts_idx_tbi_name(fp_->idx, rid,
                bcf_hdr_id2name(header_, rid))) < 0)
            {
                std::cerr << "[VcfWriter::write] Failed to get tabix tid! "
                    << std::endl;
                std::exit(1);
            }

            if (hts_idx_push(fp_->idx, tid, start, end,
                bgzf_tell(fp_->fp.bgzf), 1) < 0)
            {
                std::cerr << "[VcfWriter::write] Failed to push hts index!"
                    << std::endl;
                std::exit(1);
            }
        }
    } else {
        kputsn(v_str.c_str(), v_str.size()-1, ks_clear(&tmp_var_ks_));
        if (v_ == nullptr) v_ = bcf_init1();
        vcf_parse(&tmp_var_ks_, header_, v_);
        write(v_);
    }
}


void VcfWriter::write(bcf1_t *v)
{
    if (bcf_write(fp_, header_, v) < 0) {
        std::cerr << "[VcfWriter::write] Failed to write vcf record. "
            << v->rid << ":" << v->pos
            << std::endl;
        std::exit(1);
    }
}



