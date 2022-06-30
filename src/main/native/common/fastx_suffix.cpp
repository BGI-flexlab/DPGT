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
