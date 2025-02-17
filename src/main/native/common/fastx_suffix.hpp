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
