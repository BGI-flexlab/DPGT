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
#include "check_jemalloc.hpp"
#include "spdlog/spdlog.h"

bool check_jemalloc() {
    if (mallctl != nullptr) {
        const char *j;
        size_t s = sizeof(j);
        mallctl("version", &j,  &s, NULL, 0);
        spdlog::info("jemalloc {} is in use.", j);
        return true;
    } else {
        spdlog::warn("jemalloc not in use!");
        return false;
    }
}
