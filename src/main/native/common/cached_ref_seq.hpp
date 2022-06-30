#ifndef DPGT_CACHED_REF_SEQ_HPP
#define DPGT_CACHED_REF_SEQ_HPP

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <iostream>
#include "boost/range/detail/common.hpp"
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "sequence_dictionary.hpp"


// if we want to print debug information about cache hit efficency ?
#define DPGT_CACHED_REF_SEQ_PRINT_EFFICENCY 0


/**
 * A cached reference sequence class
 *
 * Not thread safe!
 *
 * Automatically upper-cases the bases coming in, unless the flag preserve_case
 * is explicitly set.
 * Automatically converts IUPAC bases to Ns, unless the flag preserve_IUPAC
 * is explicitly set.
 */
class CachedRefSeq {
private:
    struct Cache {
        int32_t tid = -1;
        int64_t start = -1;
        int64_t stop = -1;
        std::string seq = "";
    };

public:

    // if If we are printing efficiency info, what frequency should we do it at?
    static const bool PRINT_FREQUENCY = 10000;

    // The default cache size in bp
    static const int64_t DEFAULT_CACHE_SIZE = 1000000;

private:
    const int64_t cache_size_;
    const int64_t cache_miss_backup_;
    const bool preserve_case_;
    const bool preserve_IUPAC_;

    faidx_t *fai_ = nullptr;
    SequenceDictionary *seq_dict_ = nullptr;

    // for print cache efficency
    int cache_hits_ = 0;
    int cache_misses_ = 0;

    Cache cache_;

    int64_t calculateCacheMissBackup() const {
        return cache_size_ / 1000 > 1 ? cache_size_ / 1000 : 1;
    }

    void loadFai(const std::string &path) {
        fai_ = fai_load(path.c_str());
        if (fai_ == nullptr) {
            std::cerr << "Error! Fail to load fai index for fasta " << path
                << std::endl;
            std::exit(1);
        }
    }

    void loadSeqDict(const std::string &path) {
        seq_dict_ = new SequenceDictionary();
        seq_dict_->Load(path);
    }

    std::string fetchSeqAt(const std::string &target,
        int64_t start, int64_t stop);
    
    void transformSeq(std::string &seq);

public:
    CachedRefSeq(const std::string &path,
        int64_t cache_size = CachedRefSeq::DEFAULT_CACHE_SIZE,
        bool preserve_case = false, bool preserve_IUPAC = false);
    
    CachedRefSeq(const CachedRefSeq &other) = delete;
    CachedRefSeq(CachedRefSeq &&other) noexcept = delete;
    CachedRefSeq &operator=(const CachedRefSeq &other) = delete;
    CachedRefSeq &operator=(CachedRefSeq &&other) noexcept = delete;

    ~CachedRefSeq();
    
    /**
     * @brief print cache efficency. 100 * cache_hit_ / query_times
     */
    void printEfficency() {
        fprintf(stderr, "CacheRefSeq cache efficency: %.6f",
            100.0 * cache_hits_ / (cache_hits_ + cache_misses_));
    }

    /**
     * @brief Get complete sequence of input contig
     * 
     * @param target contig name
     */
    std::string getSequence(const std::string &target);

    /**
     * @brief Get complete sequence of input contig
     * 
     * @param tid contig index
     */
    std::string getSequence(int32_t tid);

    /**
     * @brief Get subsequence at given interval
     * 
     * @param target contig name
     * @param start 0-based start(include) of the interval
     * @param stop 0-based stop(include) of the interval
     * @return std::string 
     */
    std::string getSubsequenceAt(const std::string &target,
        int64_t start, int64_t stop);
    
    /**
     * @brief Get subsequence at given interval
     * 
     * @param tid contig index
     * @param start 0-based start(include) of the interval
     * @param stop 0-based stop(include) of the interval
     * @return std::string 
     */
    std::string getSubsequenceAt(int32_t tid, int64_t start, int64_t stop);

    int getCacheSize() const {
        return cache_size_;
    }

    int getCacheMissBackup() const {
        return cache_miss_backup_;
    }

    bool isPreserveCase() const {
        return preserve_case_;
    }

    bool isPreserveIUPAC() const {
        return preserve_IUPAC_;
    }

    faidx_t *getFai() const {
        return fai_;
    }

    SequenceDictionary *getSequenceDict() const {
        return seq_dict_;
    }
};

#endif  // DPGT_CACHED_REF_SEQ_HPP
