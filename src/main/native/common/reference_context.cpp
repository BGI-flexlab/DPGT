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
#include "common/reference_context.hpp"
#include "common/simple_interval.hpp"
#include <cstdint>



void ReferenceContext::setWindow(int64_t leading_bases, int64_t trailing_bases)
{
    if (leading_bases < 0) {
        std::cerr << "[ReferenceContext::setWindow] Error! leading bases "
            << "can not be negative" << std::endl;
        std::exit(1);
    }

    if (trailing_bases < 0) {
        std::cerr << "[ReferenceContext::setWindow] Error! trailing_bases bases"
            << " can not be negative" << std::endl;
        std::exit(1);
    }

    if (interval_.IsNull() || (leading_bases == 0 && trailing_bases == 0)) {
        window_ = interval_;
    } else {
        int64_t start = interval_.start - leading_bases;
        start = trimToContigStart(start);
        int64_t end = interval_.end + trailing_bases;
        end = trimToContigEnd(end);
        window_ = SimpleInterval(interval_.tid, start,
            end);
    }

    sequence_ = "";
}
