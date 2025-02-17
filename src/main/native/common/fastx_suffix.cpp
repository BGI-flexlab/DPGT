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
//
// Created by yangqi735 on 2021/4/22.
//

#include "common/fastx_suffix.hpp"
#include "boost/algorithm/string.hpp"



const std::vector<std::string> FastxSuffix::kFastaSuffix = {
        ".fasta", ".fa", ".fasta.gz", ".fa.gz"};

const std::vector<std::string> FastxSuffix::kFastqSuffix = {
        ".fastq", ".fq", ".fastq.gz", ".fq.gz"};

std::string FastxSuffix::GetSuffix(
        const std::string &input,
        const std::vector<std::string> &suffix_list)
{
    for (auto const &suffix: suffix_list)
    {
        if (boost::iends_with(input, suffix))
            return suffix;
    }
    return "";
}

bool FastxSuffix::IsFasta(const std::string &input) {
    return !GetSuffix(input, kFastaSuffix).empty();
}

bool FastxSuffix::IsFastq(const std::string &input) {
    return !GetSuffix(input, kFastqSuffix).empty();
}

std::string FastxSuffix::TrimSuffix(const std::string &fastx) {
    std::string suffix;
    if (!(suffix = GetSuffix(fastx, kFastaSuffix)).empty())
    {
        return fastx.substr(0, fastx.size() - suffix.size());
    } else {
        suffix = GetSuffix(fastx, kFastqSuffix);
        return fastx.substr(0, fastx.size() - suffix.size());
    }
}

std::string FastxSuffix::TrimFastaSuffix(const std::string &fasta) {
    std::string suffix = GetSuffix(fasta, kFastaSuffix);
    return fasta.substr(0, fasta.size() - suffix.size());
}

std::string FastxSuffix::TrimFastqSuffix(const std::string &fastq) {
    std::string suffix = GetSuffix(fastq, kFastaSuffix);
    return fastq.substr(0, fastq.size() - suffix.size());
}
