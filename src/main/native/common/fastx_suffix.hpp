//
// Created by yangqi735 on 2021/4/22.
//

#ifndef FASTX_SUFFIX_HPP
#define FASTX_SUFFIX_HPP

#include <vector>
#include <string>


/**
 * Used for operate suffix of fasta or fastq file.
 */
class FastxSuffix {
public:
    static const std::vector<std::string> kFastaSuffix;
    static const std::vector<std::string> kFastqSuffix;

    static std::string GetSuffix(
            const std::string &input,
            const std::vector<std::string> &suffix_list);

    static bool IsFasta(const std::string &input);
    static bool IsFastq(const std::string &input);


    static std::string TrimSuffix(const std::string &fastx);
    static std::string TrimFastaSuffix(const std::string &fasta);
    static std::string TrimFastqSuffix(const std::string &fastq);
};


#endif //FASTX_SUFFIX_HPP
