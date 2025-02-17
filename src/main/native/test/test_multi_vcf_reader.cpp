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

