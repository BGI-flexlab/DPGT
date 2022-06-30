#include <cstddef>
#include <cstdint>
#include "vcf/vcf_reader.hpp"
#include "variant_site_finder.hpp"


VariantSiteSet FindVariantSite(const std::vector<std::string> &vcfpaths,
    const char *chrom, int64_t start, int64_t end)
{
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
                record_end = std::min(end, record_end);
                for (; j < record_end; ++j) vs_set.set(j);
            }
        }
        bcf_destroy(record);
    }
    return vs_set;
}
