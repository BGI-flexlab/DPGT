#include "common/fastx_seq.hpp"


std::string FastxSeq::Str() const
{
    std::ostringstream ostr;
    if (qual.empty())  // fasta
    {
        ostr << ">" << name;
        if (!comment.empty())
            ostr << " " << comment;
        ostr << "\n" << seq << std::endl;
    } else  // fastq
    {
        ostr << "@" << name;
        if (!comment.empty())
            ostr << " " << comment;
        ostr << "\n" << seq << "\n+\n" << qual << std::endl;
    }
    return ostr.str();
}


int FastxReader::Read(FastxSeq &fastx_seq)
{
    int r;
    if ( (r = kseq_read(seq)) >= 0 )
    {
        fastx_seq = FastxSeq(seq);
        return r;
    } else if (r == -1)
    {
        return r; // reach to the end-of-file
    } else if (r == -2)
    {
        fprintf(stderr, "[FastxSeq::%s] Error! Truncated quality string!\n",
            __func__);
        std::exit(-2);
    } else
    {
        fprintf(stderr, "[FastxSeq::%s] Error! Error reading stream!\n",
            __func__);
        std::exit(-3);
    }
}


int FastxReader::ReadUntil(std::vector<FastxSeq> &fastx_seqs, int n)
{
    int r = 0;
    int i = 0;
    while (i < n && (r = kseq_read(seq)) >= 0)
    {
        fastx_seqs.push_back(FastxSeq(seq));
        ++i;
    }

    if (r == -2)
    {
        fprintf(stderr, "[FastxSeq::%s] Error! Truncated quality string!\n",
            __func__);
        std::exit(-2);
    } 

    if (r == -3)
    {
        fprintf(stderr, "[FastxSeq::%s] Error! Error reading stream!\n",
            __func__);
        std::exit(-3);
    }

    return i;
}

