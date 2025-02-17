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
#include <cstring>
#include <iostream>
#include <sstream>
#include "utils.hpp"



void Utils::validateArg(bool condition, const std::string &msg) {
    if (!condition) {
        std::cerr << msg << std::endl;
        std::exit(1);
    }
}


std::string Utils::repeatJoin(
    const std::string &u, const std::string &s, int n)
{
    std::ostringstream out_str;
    int last = n - 1;
    for (int i = 0; i < last; ++i) {
        out_str << u << s;
    }
    out_str << u;
    return out_str.str();
}


void Utils::krepeatPuts(kstring_t *o, const std::string &u,
    const std::string &s, int n)
{
    int last = n - 1;
    for (int i = 0; i < last; ++i) {
        kputs(u.c_str(), o);
        kputs(s.c_str(), o);
    }
    kputs(u.c_str(), o);
}
