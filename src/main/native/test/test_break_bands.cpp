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
#include <cassert>
#include <cstddef>
#include <iostream>
#include <set>
#include <cstdint>
#include <vector>

std::set<int64_t> getIntermediateStopSites(
    int64_t start, int64_t end, int64_t break_bands_multiple)
{
    std::set<int64_t> sites_to_stop;
    if (break_bands_multiple > 0) {
        int64_t block_end_pos = 
            ((start + 1) / break_bands_multiple + 1) *
            break_bands_multiple - 1;
        for (;block_end_pos <= end;
            block_end_pos += break_bands_multiple)
        {
            sites_to_stop.insert(block_end_pos);
        }
    }
    return sites_to_stop;
}

int main() {
    int64_t start = 25;
    int64_t end = 48;
    int64_t break_bands_multiple = 4;
    std::set<int64_t> sites_to_stop = getIntermediateStopSites(start, end, break_bands_multiple);
    std::vector<int64_t> expected = {27, 31, 35, 39, 43, 47};
    assert(sites_to_stop.size() == expected.size());
    int i = 0;
    for (auto s: sites_to_stop) {
        assert(s == expected[i++]);
        std::cout << s << std::endl;
    }

    std::cout << "####" << std::endl;

    start = 0;
    end = 5;
    break_bands_multiple = 1;
    sites_to_stop = getIntermediateStopSites(start, end, break_bands_multiple);
    expected = {1, 2, 3, 4, 5};
    assert(sites_to_stop.size() == expected.size());
    i = 0;
    for (auto s: sites_to_stop) {
        assert(s == expected[i++]);
        std::cout << s << std::endl;
    }
}
