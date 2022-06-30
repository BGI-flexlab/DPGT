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
