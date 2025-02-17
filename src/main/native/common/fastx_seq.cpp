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

