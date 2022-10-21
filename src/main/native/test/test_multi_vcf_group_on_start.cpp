#include <cstdio>
#include <vector>
#include <limits>

#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "vcf/multi_vcf_reader.hpp"
#include "vcf/variant_context.hpp"


static void PrintVariants(const std::vector<VariantContext> records,
    const bcf_hdr_t *h)
{
    kstring_t line = {0, 0, NULL};
    for (auto r: records) {
        vcf_format1(h, r.getVariant(), ks_clear(&line));
        printf("%s", line.s);
    }
    if (line.s) free(line.s);
}

static bool CheckGroup(const std::vector<VariantContext> records) {
    int rid = records.front().getContig();
    int pos = records.front().getStart();
    for (auto r: records) {
        if (r.getContig() != rid || r.getStart() != pos)
        {
            std::cerr << "Error! False grouped variants! " << rid << ":" << pos << "\t" << r.getContig() << ":" << r.getStart() << std::endl;
            return false;
        }
    }
    return true;
}


class MultiVariantGroupOnStart {
public:
    std::vector<VariantContext> current_variants;
    int first_current_var_start = 0;
    int last_current_var_start = 0;

    MultiVcfReader *reader;

    explicit MultiVariantGroupOnStart(MultiVcfReader *reader_): reader(reader_) {}

    ~MultiVariantGroupOnStart() {}

    /**
     * @brief 输入1个新变异记录，如果新变异记录的开始位置与已经分组的变异
     * 记录不同，那么打印保存的变异记录
     * 
     */
    void apply(const VariantContext &variant_context);

    void run() {
        VariantContext variant_context;
        while ( !(variant_context = reader->Read()).isNull() ) {
            apply(variant_context);
        }
        if (!current_variants.empty()) {
            CheckGroup(current_variants);
            PrintVariants(current_variants, reader->header());
        }
    }
};


void MultiVariantGroupOnStart::apply(const VariantContext &variant_context) {
    if (current_variants.empty()) {
        first_current_var_start = variant_context.getStart();
    } else if (current_variants.front().getContig() != variant_context.getContig() ||
        last_current_var_start < variant_context.getStart() - 0 ||
        first_current_var_start < -std::numeric_limits<int>::max())
    {
        CheckGroup(current_variants);
        PrintVariants(current_variants, reader->header());
        printf("\n====================\n\n");
        // for (auto r: current_variants) bcf_destroy(r);
        current_variants.clear();
        first_current_var_start = variant_context.getStart();
    }

    current_variants.push_back(variant_context);
    last_current_var_start = variant_context.getStart();
}


int main(int argc, char **argv) {
    MultiVcfReader reader({argv[1], argv[2]}, false);

    MultiVariantGroupOnStart grouper(&reader);

    kstring_t header_str = {0, 0, NULL};
    bcf_hdr_format(grouper.reader->header(), false, &header_str);
    printf("%s\n", header_str.s);

    grouper.run();
}



