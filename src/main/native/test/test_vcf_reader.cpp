#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "common/simple_interval.hpp"
#include "htslib/vcf.h"
#include "htslib/kstring.h"
#include "vcf/vcf_reader.hpp"
#include "vcf/vcf_ibuffer.hpp"


int main(int argc, char **argv) {
    VcfReader reader(argv[1], true);

    VcfIBuffer buffer(&reader);

    std::vector<SimpleInterval> intervals{SimpleInterval(19, 10436226, 10436247)};
    reader.QueryIntervals(intervals);

    bcf1_t *record;
    kstring_t line = {0, 0, NULL};
    while ( (record = buffer.Read()) ) {
        vcf_format1(reader.header(), record, ks_clear(&line));
        printf("%s", line.s);
    }

    if (line.s) free(line.s);
}

