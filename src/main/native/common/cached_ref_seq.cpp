#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <algorithm>
#include "cached_ref_seq.hpp"
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "vcf/nucleotide.hpp"



CachedRefSeq::CachedRefSeq(const std::string &path, int64_t cache_size,
    bool preserve_case, bool preserve_IUPAC):
    cache_size_(cache_size),
    cache_miss_backup_(calculateCacheMissBackup()),
    preserve_case_(preserve_case), preserve_IUPAC_(preserve_IUPAC)
{
    loadFai(path);
    loadSeqDict(path);
}

CachedRefSeq::~CachedRefSeq() {
    if (fai_) fai_destroy(fai_);
    if (seq_dict_) delete seq_dict_;
}

std::string CachedRefSeq::fetchSeqAt(const std::string &target,
    int64_t start, int64_t stop)
{
    int64_t target_len = 0;
    char *seq = faidx_fetch_seq64(fai_, target.c_str(), start,
        stop, &target_len);
    if (target_len < 0 || seq == nullptr) {
        std::cerr << "[CachedRefSeq::fetchSeqAt] Error! Can not fetch "
                        "sequence! target=" << target.c_str() << std::endl;
        free(seq);
        std::exit(1);
    }
    std::string res = seq;  // copy seq to c++ string
    free(seq);
    return res;
}

void CachedRefSeq::transformSeq(std::string &seq) {
    if (!preserve_case_ && !preserve_IUPAC_) {
        std::transform(seq.begin(), seq.end(),
            seq.begin(), nucleotide::ToUpperACTGN);
    } else if (!preserve_case_) {
        std::transform(seq.begin(), seq.end(),
            seq.begin(), toupper);
    } else if (!preserve_IUPAC_) {
        auto iupac_to_N = [](char c) -> char {
            if (nucleotide::IsStandardNt(c)) return c;
            return 'N';
        };
        std::transform(seq.begin(), seq.end(),
            seq.begin(), iupac_to_N);
    }
}

std::string CachedRefSeq::getSequence(const std::string &target) {
    return fetchSeqAt(target, 0, HTS_POS_MAX);
}

std::string CachedRefSeq::getSequence(int32_t tid) {
    std::string target;
    if (seq_dict_->SequenceIndexToName(tid, target) != 0) {
        std::cerr << "[CachedRefSeq::getSequence] Error! Can not convert"
            " tid to target name. tid=" << tid << std::endl;
        std::exit(1);
    }
    return getSequence(target);
}

std::string CachedRefSeq::getSubsequenceAt(
    const std::string &target, int64_t start, int64_t stop)
{
    std::string result;

    int32_t tid;
    if (seq_dict_->SequenceNameToIndex(target, tid) != 0) {
        std::cerr << "[CachedRefSeq::getSequence] Error! Can not convert"
            " target name to tid. target=" << tid << std::endl;
        std::exit(1);
    }

    if (stop - start + 1 > cache_size_) {
        ++cache_misses_;
        result = fetchSeqAt(target, start, stop);
        transformSeq(result);
    } else {

        int64_t target_len;
        seq_dict_->SequenceIndexToLength(tid, target_len);

        if (start < cache_.start || stop > cache_.stop ||
            cache_.seq.empty() || cache_.tid != tid)
        {
            ++cache_misses_;
            cache_.tid = tid;
            cache_.start = std::max(start - cache_miss_backup_, 0L);
            cache_.stop = std::min(start + cache_size_ + cache_miss_backup_,
                target_len-1);
            cache_.seq = fetchSeqAt(target, cache_.start, cache_.stop);

            transformSeq(cache_.seq);
        } else {
            ++cache_hits_;
        }

        result = cache_.seq.substr(start - cache_.start, stop - start + 1);
    }

#if DPGT_CACHED_REF_SEQ_PRINT_EFFICENCY
    if ((cache_hits_ + cache_misses_) % PRINT_FREQUENCY == 0)
    {
        printEfficency();
    }
#endif // DPGT_CACHED_REF_SEQ_PRINT_EFFICENCY

    return result;
}


std::string CachedRefSeq::getSubsequenceAt(int32_t tid,
    int64_t start, int64_t stop)
{
    std::string target;
    if (seq_dict_->SequenceIndexToName(tid, target) != 0) {
        std::cerr << "[CachedRefSeq::getSubsequenceAt] Error! Can not convert"
            " tid to target name. tid=" << tid << std::endl;
        std::exit(1);
    }
    return getSubsequenceAt(target, start, stop);
}
