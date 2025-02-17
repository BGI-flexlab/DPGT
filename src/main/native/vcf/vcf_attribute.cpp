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
#include "vcf/vcf_attribute.hpp"



const std::map<int, std::string> VcfAttributeBase::TYPE_STRINGS = {
    {BCF_HT_FLAG, "FLAG"},
    {BCF_HT_INT, "INT"},
    {BCF_HT_REAL, "FLOAT"},
    {BCF_HT_STR, "STRING"},
    {BCF_HT_LONG, "LONG"}
};


const std::array<std::string, 5> VcfAttributeBase::SIZE_TYPE_STRINGS = {
    "FIXED", "VAR", "A", "G", "R"};
