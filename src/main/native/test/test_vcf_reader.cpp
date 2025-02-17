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
#include <vector>
#include "common/simple_interval.hpp"
#include "htslib/vcf.h"
#include "htslib/kstring.h"
#include "vcf/vcf_reader.hpp"
#include "vcf/vcf_ibuffer.hpp"


int main(int argc, char **argv) {
    VcfReader reader(argv[1], true);

    // reader.Querys("chr20", 10436226, 10436247);

    VcfIBuffer buffer(&reader);

    std::vector<SimpleInterval> intervals{SimpleInterval(19, 10436226, 10436247)};
    buffer.QueryIntervals(intervals);

    bcf1_t *record;
    kstring_t line = {0, 0, NULL};
    while ( (record = buffer.Read()) ) {
        vcf_format1(reader.header(), record, ks_clear(&line));
        printf("%s", line.s);
    }

    if (line.s) free(line.s);
}

