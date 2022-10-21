#include <cstdio>
#include <iostream>
#include <cstdlib>
#include "htslib/vcf.h"
#include "htslib/kstring.h"
#include "vcf/vcf_reader.hpp"
#include "vcf/multi_vcf_reader.hpp"


int main(int argc, char **argv) {
    MultiVcfReader reader({argv[1], argv[2]}, false);

    // reader.Querys("chr1", 0, 14522);

    bcf1_t *record;
    bcf_hdr_t *header = reader.header();
    kstring_t line = {0, 0, NULL};
    while ( (record = reader.Read()) ) {
        vcf_format1(header, record, ks_clear(&line));
        printf("%s", line.s);
    }

    bcf_hdr_destroy(header);

    if (line.s) free(line.s);
}

