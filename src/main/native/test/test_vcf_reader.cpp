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
    // VcfReader reader(argv[1], argv[2], true, false);
    VcfReader reader(argv[1], argv[2], true, true);

    // reader.Querys("chr20", 10436226, 10436247);

    // VcfIBuffer buffer(&reader);
    std::vector<SimpleInterval> intervals{
        SimpleInterval(19, 1000000, 2000000),
        SimpleInterval(19, 3000000, 3500000)};
    // buffer.QueryIntervals(intervals);
    reader.QueryIntervals(intervals);
    bcf1_t *record = bcf_init1();
    kstring_t line = {0, 0, NULL};
    while ( reader.Read(record) ) {
        vcf_format1(reader.header(), record, ks_clear(&line));
        printf("%s", line.s);
    }

    bcf_destroy1(record);
    if (line.s) free(line.s);
}

