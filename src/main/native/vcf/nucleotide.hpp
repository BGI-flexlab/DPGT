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
#ifndef DPGT_NUCLEOTIDE_HPP
#define DPGT_NUCLEOTIDE_HPP

#include <string>
#include "htslib/hts.h"

namespace nucleotide {

inline
bool IsEqual(char n1, char n2)
{
    return seq_nt16_table[(int)n1] == seq_nt16_table[(int)n2];
}

/**
 * Is the char a standard nucleotide(in ACGTacgt) ?
 * @param n input char
 * @return true if the char is a standard nucleotide, otherwise return false.
 */
inline
bool IsStandardNt(char n) {
    // "=ACMGRSVTWYHKDBN";
    static const bool is_standard_nt16[] =
            {0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0};
           //=  A  C  M  G  R  S  V  T  W  Y  H  K  D  B  N
    return is_standard_nt16[seq_nt16_table[(int)n]];
}


/**
 * Is the char a standard nucleotide include 'Nn' (in ACGTNacgtn) ?
 * @param n input char
 * @return true if the char is a standard nucleotide, otherwise return false.
 */
inline
bool IsStandardNtIncludeNn(char n) {
    // "=ACMGRSVTWYHKDBN";
    static const bool is_standard_nt16_nn[] =
            {0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
           //=  A  C  M  G  R  S  V  T  W  Y  H  K  D  B  N
    return is_standard_nt16_nn[seq_nt16_table[(int)n]];
}

/**
 * 1. transform base to upper case
 * 2. transform non-ACGTN IUPAC ambiguity code to 'N' 
 */
inline
char ToUpperACTGN(char n) {
    static const char seq_nt16_str_abn[] = "NACNGNNNTNNNNNNN";
    return seq_nt16_str_abn[seq_nt16_table[(int)n]];
}

inline
std::string ToStr(char n1)
{
    return std::string({seq_nt16_str[seq_nt16_table[(int)n1]]});
}

}  // nucleotide



#endif  // DPGT_NUCLEOTIDE_HPP
