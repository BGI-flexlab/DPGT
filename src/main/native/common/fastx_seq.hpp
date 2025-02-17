/**
This file is part of DPGT.
Copyright (C) 2022 BGI.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// License End
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
