#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "common/simple_interval.hpp"
#include "htslib/vcf.h"
#include "zlib.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "vcf/vcf_writer.hpp"
#include "vcf/vcf_reader.hpp"


int main(int argc, char **argv) {
    BGZF *infp = bgzf_open(argv[1], "r");

    VcfReader reader{argv[1]};
    bcf_hdr_t *header = reader.header();

    VcfWriter writer{argv[2]};
    writer.writeHeader(header);

    kstring_t line = {0, 0, NULL};
    int r;
    while ((r = bgzf_getline(infp, '\n', &line) >= 0)) {
        int j = 0, n = 0, pos = 0, rlen = 0;
        if (line.s[0] == '#') continue;
        for (size_t i = 0; i < line.l; ++i)
        {
            if (line.s[i] == '\t') {
                n += 1;
                if (n == 2) {
                    // POS
                    line.s[i] = '\0';
                    pos = atol(line.s+j) - 1;
                    line.s[i] = '\t';
                }
                if (n == 4) {
                    // REF
                    rlen = i - j;
                    break;
                }
                j = i + 1;
            }
        }
        kputc('\n', &line);
        writer.write(ks_c_str(&line), 19, pos, pos + rlen);
    }

    bgzf_close(infp);
}