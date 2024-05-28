#ifndef DPGT_VARIANT_SITE_FINDER_HPP
#define DPGT_VARIANT_SITE_FINDER_HPP

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include "boost/dynamic_bitset/dynamic_bitset.hpp"
#include "variant/variant_site_set.hpp"


VariantSiteSet FindVariantSite(const std::vector<std::string> &vcfpaths,
    const std::vector<std::string> &indexpaths, bool use_lix,
    const char *chrom, int64_t start, int64_t end);


#endif  // DPGT_VARIANT_SITE_FINDER_HPP
