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
#include <cstddef>
#include <cstdint>
#include "vcf/vcf_reader.hpp"
#include "variant_site_finder.hpp"


VariantSiteSet FindVariantSite(const std::vector<std::string> &vcfpaths,
    const char *chrom, int64_t start, int64_t end)
{
    int64_t end1 = end + 1;
    VariantSiteSet vs_set{start, end};
    for (size_t i = 0; i < vcfpaths.size(); ++i) {
        VcfReader reader{vcfpaths[i], true};
        reader.Querys(chrom, start, end);
        bcf1_t *record = bcf_init1();
        while (reader.Read(record)) {
            // GATK gvcf record has two alleles(ref, <NON-REF>) for non-variant site,
            // if number of alleles greater than 2, then there must be at least one
            // true alt allele
            if (record->n_allele > 2) {
                int64_t record_end = record->pos + record->rlen;  // end(0-based) + 1
                int64_t j = std::max(start, record->pos);
                record_end = std::min(end1, record_end);
                for (; j < record_end; ++j) vs_set.set(j);
            }
        }
        bcf_destroy(record);
    }
    return vs_set;
}
