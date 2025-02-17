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
#ifndef DPGT_CHECK_JEMALLOC_HPP
#define DPGT_CHECK_JEMALLOC_HPP

#include <stddef.h>

extern "C"
{
    // weak symbol: resolved at runtime by the linker if we are using jemalloc, nullptr otherwise
    int mallctl(const char *name, void *oldp, size_t *oldlenp, void *newp, size_t newlen) __attribute__((weak));
}


bool check_jemalloc();


#endif  // DPGT_CHECK_JEMALLOC_HPP
