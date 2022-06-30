#ifndef FASTX_SEQ_HPP
#define FASTX_SEQ_HPP


#include <string>
#include <sstream>
#include <cstring>
#include <vector>

#include "zlib.h"
#include "htslib/kseq.h"

KSEQ_INIT(gzFile, gzread)


/**
 * @brief fasta/fastq record class
 */
struct FastxSeq
{
    FastxSeq() = default;
    FastxSeq(const kseq_t *seq_)
    {
        // fprintf(stdout, "%s\t%s\n", seq_->name.s, seq_->seq.s);
        name = seq_->name.s;
        if (seq_->comment.l > 0) comment = seq_->comment.s;
        seq = seq_->seq.s;
        if (seq_->qual.l > 0) qual = seq_->qual.s;
    }
    ~FastxSeq() {};

    std::string Str() const;
    
    std::string name = "";
    std::string comment = "";
    std::string seq = "";
    std::string qual = "";
};


/**
 * @brief a fasta/fastq reader using kseq
 */
class FastxReader
{
public:
    FastxReader() = delete;
    FastxReader(const std::string &fn)
    {
        fp = gzopen(fn.c_str(), "r");
        seq = kseq_init(fp);
    }

    // disable copy and move construct/assignment
    FastxReader(const FastxReader &other) = delete;
    FastxReader(FastxReader &&other) = delete;
    FastxReader &operator=(const FastxReader &other) = delete;
    FastxReader &operator=(FastxReader &&other) = delete;

    ~FastxReader()
    {
        kseq_destroy(seq);
        gzclose(fp);
    }

    // read 1 sequence
    int Read(FastxSeq &fastx_seq);
    
    // read n sequence or read until the end of file
    int ReadUntil(std::vector<FastxSeq> &fastx_seqs, int n);

private:
    gzFile fp;
    kseq_t *seq;
};


#endif  // FASTX_SEQ_HPP
