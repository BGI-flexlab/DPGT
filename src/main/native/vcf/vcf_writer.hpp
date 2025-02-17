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
#ifndef DPGT_VCF_WRITER_HPP
#define DPGT_VCF_WRITER_HPP

#include <cstdint>
#include <iostream>
#include <string>
#include "htslib/vcf.h"


/**
 * @brief vcf/bcf wruter based on htslib
 */
class VcfWriter {
private:
    htsFile *fp_ = nullptr;
    std::string idx_fn_ = "";
    bcf_hdr_t *header_ = nullptr;

    bcf1_t *v_ = nullptr;
    kstring_t tmp_var_ks_ = {0, 0, NULL};

public:
    VcfWriter() = default;
    VcfWriter(const std::string &fn);
    VcfWriter(const VcfWriter &other) = delete;
    VcfWriter(VcfWriter &&other) noexcept = delete;
    VcfWriter &operator=(const VcfWriter &other) = delete;
    VcfWriter &operator=(VcfWriter &&other) noexcept = delete;
    ~VcfWriter() { close(); }

    void open(const std::string &fn);

    void close();

    void writeHeader(bcf_hdr_t *hdr);

    void write(kstring_t *v_str,
        int32_t rid, int64_t start, int64_t end);
    void write(bcf1_t *v);
};


#endif  // DPGT_VCF_WRITER_HPP
