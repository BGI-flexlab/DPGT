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
#ifndef DPGT_UTILS_HPP
#define DPGT_UTILS_HPP


#include <string>
#include "htslib/kstring.h"

class Utils {
public:
    static void validateArg(bool condition, const std::string &msg);

    static
    std::string repeatJoin(const std::string &u, const std::string &s, int n);

    static void krepeatPuts(kstring_t *o, const std::string &u,
        const std::string &s, int n);
};


template <typename PtrT>
struct PtrLess {
    bool operator()(PtrT lhs, PtrT rhs) const
    {
        return *lhs < *rhs;
    }
};


template <typename T>
struct ArrayBuff {
    int l; // length
    int m; // capacity
    T *data;
};

#endif  // DPGT_UTILS_HPP
