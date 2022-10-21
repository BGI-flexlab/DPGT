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
