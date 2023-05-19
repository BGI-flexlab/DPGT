#include "vcf/vcf_reader.hpp"
#include "common/simple_interval.hpp"
#include "htslib/tbx.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include <cstdint>
#include <cstdlib>
#include <iostream>


#define MAX_CSI_COOR ((1LL << (14 + 30)) - 1)


VcfReader::VcfReader(VcfReader &&reader) noexcept
{
    // move by copy pointers and values
    this->file_ptr_ = reader.file_ptr_;
    this->tbx_idx_ = reader.tbx_idx_;
    this->bcf_idx_ = reader.bcf_idx_;
    this->header_ = reader.header_;
    this->itr_ = reader.itr_;
    this->streaming_ = reader.streaming_;
    this->tmps_ = reader.tmps_;

    // set other reader to null
    reader.file_ptr_ = nullptr;
    reader.tbx_idx_ = nullptr;
    reader.bcf_idx_ = nullptr;
    reader.header_ = nullptr;
    reader.itr_ = nullptr;
    reader.tmps_ = {0, 0, NULL};
}


VcfReader &VcfReader::operator==(VcfReader &&reader) noexcept {
    if (&reader != this) {
        // free resources of this
        Close();

        // move by copy pointers and values
        this->file_ptr_ = reader.file_ptr_;
        this->tbx_idx_ = reader.tbx_idx_;
        this->bcf_idx_ = reader.bcf_idx_;
        this->header_ = reader.header_;
        this->itr_ = reader.itr_;
        this->streaming_ = reader.streaming_;
        this->tmps_ = reader.tmps_;

        // set other reader to null
        reader.file_ptr_ = nullptr;
        reader.tbx_idx_ = nullptr;
        reader.bcf_idx_ = nullptr;
        reader.header_ = nullptr;
        reader.itr_ = nullptr;
        reader.tmps_ = {0, 0, NULL};
    }
    return *this;
}


void VcfReader::Open(const std::string &file_name, bool require_index) {
    char fmode[5];
    strcpy(fmode, "r");
    vcf_open_mode(fmode+1, file_name.c_str(), NULL);
    file_ptr_ = hts_open(file_name.c_str(), fmode);
    if ( ! file_ptr_ ) {
        std::cerr << "[VcfReader::Open] Error! Fail to open " << file_name
            << std::endl;
        std::exit(1);
    }

    if (require_index) {
        streaming_ = false;
        if ( file_ptr_->format.format==vcf )
        {
            if ( file_ptr_->format.compression!=bgzf )
            {
                std::cerr << "[VcfReader::Open] Error! Require index but"
                    << " input file is not a bgzip vcf. " << file_name
                    << std::endl;
                std::exit(1);
            }

            tbx_idx_ = tbx_index_load(file_name.c_str());
            if ( !tbx_idx_ )
            {
                std::cerr << "[VcfReader::Open] Error! Failed to load tbx "
                    << "index. " << file_name << std::endl;
                std::exit(1);
            }

            header_ = bcf_hdr_read(file_ptr_);
        }
        else if ( file_ptr_->format.format==bcf )
        {
            if ( file_ptr_->format.compression!=bgzf )
            {
                std::cerr << "[VcfReader::Open] Error! Require index but"
                    << " input file is not a bgzip bcf. " << file_name
                    << std::endl;
                std::exit(1);
            }

            bcf_idx_ = bcf_index_load(file_name.c_str());
            if ( !bcf_idx_ )
            {
                std::cerr << "[VcfReader::Open] Error! Failed to load bcf"
                    << " index. " << file_name << std::endl;
                std::exit(1);
            }

            header_ = bcf_hdr_read(file_ptr_);
        }
        else
        {
            std::cerr << "[VcfReader::Open] Error! Wrong input file format. "
                << file_name << std::endl;
            std::exit(1);
        }
        // seek to the start offset before any explicit call of Query*
        // Queryi(HTS_IDX_START, 0, 0);
    } else {
        if ( file_ptr_->format.format==bcf ||
            file_ptr_->format.format==vcf )
        {
            header_ = bcf_hdr_read(file_ptr_);
        }
        else
        {
            std::cerr << "[VcfReader::Open] Error! Wrong input file format. "
                << file_name << std::endl;
            std::exit(1);
        }
    }
}


VcfReader &VcfReader::Queryi(int32_t tid, int64_t start, int64_t end) {
    if (!file_ptr_) {
        std::cerr << "[VcfReader::Queryi] Error! Can not query on reader "
            << "not opened." << std::endl;
        std::exit(1);
    }

    if (streaming_) {
        std::cerr << "[VcfReader::Queryi] Error! This reader is not support"
            << " random query." << std::endl;
        std::exit(1);
    }

    if ( end>=MAX_CSI_COOR )
    {
        end = MAX_CSI_COOR;
    }
    if ( itr_ )
    {
        hts_itr_destroy(itr_);
        itr_ = nullptr;
    }
    if ( tbx_idx_ )
    {
        if (tid == HTS_IDX_START) {
            itr_ = tbx_itr_queryi(tbx_idx_, tid, start, end+1);
        } else {
            // convert tid from bcf header tid(index of chromosome in bcf
            // header) to tabix tid(index in bgzip(first column index))
            const char *name = bcf_hdr_id2name(header_, tid);
            if (name == NULL) {
                std::cerr << "[VcfReader::Queryi] Error! Failed to convert bcf"
                    " tid:" << tid << " to name" << std::endl;
                std::exit(1);
            }
            tid = tbx_name2id(tbx_idx_, name);
            if (tid < 0) {
                itr_ = nullptr;
            }
            itr_ = tbx_itr_queryi(tbx_idx_, tid, start, end+1);
        }
    }
    else
    {
        itr_ = bcf_itr_queryi(bcf_idx_, tid, start, end+1);
    }
    return *this;
}


VcfReader &VcfReader::Querys(
    const std::string &chrom, int64_t start, int64_t end)
{
    if (!file_ptr_) {
        std::cerr << "[VcfReader::Querys] Error! Can not query on reader "
            << "not opened." << std::endl;
        std::exit(1);
    }

    if (streaming_) {
        std::cerr << "[VcfReader::Querys] Error! This reader is not support"
            << " random query." << std::endl;
        std::exit(1);
    }

    int32_t tid;
    // if (tbx_idx_) {
    //     tid = tbx_name2id(tbx_idx_, chrom.c_str());
    //     if (tid == -1) {
    //         std::cerr << "[VcfReader::Querys] Error! chrom not present in "
    //             << "this file." << std::endl;
    //         std::exit(1);
    //     }
    // } else {
    tid = bcf_hdr_name2id(header_, chrom.c_str());
    if (tid == -1) {
        std::cerr << "[VcfReader::Querys] Error! chrom " << chrom
            << " not present in vcf/bcf header " << std::endl;
        std::exit(1);
    }
    // }

    return Queryi(tid, start, end);
}


VcfReader &VcfReader::QueryIntervals(
    const std::vector<SimpleInterval> &intervals)
{
    if (intervals.empty()) return *this;

    if (!file_ptr_) {
        std::cerr << "[VcfReader::QueryIntervals] Error! Can not query on "
            << "reader not opened." << std::endl;
        std::exit(1);
    }

    if (streaming_) {
        std::cerr << "[VcfReader::QueryIntervals] Error! This reader is not "
            << "support random query." << std::endl;
        std::exit(1);
    }

    intervals_ = &intervals;

    intervals_iter_ = intervals_->begin();
    const SimpleInterval &first = *intervals_iter_;
    Queryi(first.tid, first.start, first.end);
    ++intervals_iter_;

    return *this;
}


bcf1_t *VcfReader::Read(bcf1_t *record) {
    if (!file_ptr_) {
        std::cerr << "[VcfReader::Read] Error! Can not read from reader "
            << "not opened." << std::endl;
        std::exit(1);
    }

    int ret;
    if (streaming_ || itr_ == nullptr) {
        if ( file_ptr_->format.format==vcf )
        {
            if ( (ret=hts_getline(file_ptr_, KS_SEP_LINE, &tmps_)) < 0 )
                return nullptr;
            ret = vcf_parse1(&tmps_, header_, record);
            if ( ret < 0 ) {
                std::cerr << "[VcfReader::Read] Error! Failed to parse vcf line"
                    << "." << std::endl;
                std::exit(1);
            }
        }
        else if ( file_ptr_->format.format==bcf )
        {
            ret = bcf_read1(file_ptr_, header_, record);
            if ( ret < -1 ) {
                std::cerr << "[VcfReader::Read] Error! Failed to parse bcf line"
                    << "." << std::endl;
                std::exit(1);
            }
            if ( ret < 0 ) return nullptr; // no more lines or an error
        } else {
            std::cerr << "[VcfReader::Read] Error! Wrong input file format."
                << std::endl;
            std::exit(1);
        }
    } else if ( tbx_idx_ )  // bgzip vcf file
    {
        while ((ret=tbx_itr_next(file_ptr_, tbx_idx_, itr_, &tmps_)) < 0) {
            if (intervals_ != nullptr && intervals_iter_ != intervals_->end()) {
                const SimpleInterval &interval = *intervals_iter_;
                Queryi(interval.tid, interval.start, interval.end);
                ++intervals_iter_;
            } else {
                return nullptr;  // end of query interval
            }
        }
        ret = vcf_parse1(&tmps_, header_, record);
        if ( ret < 0 ) {
            std::cerr << "[VcfReader::Read] Error! Failed to parse vcf line"
                << "." << std::endl;
            std::exit(1);
        }
    } else {  // bcf file
        while ((ret = bcf_itr_next(file_ptr_, itr_, record)) < 0) {
            if ( ret < -1 ) {
                std::cerr << "[VcfReader::Read] Error! Failed to parse bcf line"
                    << "." << std::endl;
                std::exit(1);
            }
            if (intervals_ != nullptr && intervals_iter_ != intervals_->end()) {
                const SimpleInterval &interval = *intervals_iter_;
                Queryi(interval.tid, interval.start, interval.end);
                ++intervals_iter_;
            } else {
                return nullptr;  // end of query interval
            }
        }
        bcf_subset_format(header_, record);
    }
    return record;
}
