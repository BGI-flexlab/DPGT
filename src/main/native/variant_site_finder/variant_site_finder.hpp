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
#ifndef DPGT_VARIANT_SITE_FINDER_HPP
#define DPGT_VARIANT_SITE_FINDER_HPP

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include "boost/dynamic_bitset/dynamic_bitset.hpp"
#include "variant/variant_site_set.hpp"


VariantSiteSet FindVariantSite(const std::vector<std::string> &vcfpaths,
    const char *chrom, int64_t start, int64_t end);


#endif  // DPGT_VARIANT_SITE_FINDER_HPP
