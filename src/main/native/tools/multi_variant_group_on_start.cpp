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
#include "tools/multi_variant_group_on_start.hpp"
#include "common/reference_context.hpp"
#include "common/simple_interval.hpp"
#include "vcf/variant_context.hpp"
#include <cstdint>
#include <limits>


void MultiVariantGroupOnStart::apply(std::shared_ptr<VariantContext> &&var) {
    const int64_t var_start = var->getStart();

    if (ignore_intervals_outside_start_ && isWithinInterval(
        SimpleInterval(var->getContig(), var_start, var_start))) {
        return;  // skip input var
    }

    if (current_variants_.empty()) {
        first_current_var_start_ = var_start;
    } else if (var->getContig() != current_variants_.front()->getContig() ||
        last_current_var_start_ < var_start)
    {
        apply(current_variants_);
        current_variants_.clear();
        first_current_var_start_ = var_start;
    }

    current_variants_.push_back(std::move(var));
    last_current_var_start_ = var_start;
}


void MultiVariantGroupOnStart::run()
{
    std::shared_ptr<VariantContext> variant_context;
    while ((variant_context = reader_->Read()) != nullptr) {
        apply(std::move(variant_context));
    }
    finalize();
}


ReferenceContext MultiVariantGroupOnStart::makeSpanningReferenceContext(
    std::vector<std::shared_ptr<VariantContext>> &variants,
    int64_t reference_window_padding) const
{
    int32_t contig = variants.front()->getContig();
    int64_t start = std::numeric_limits<int64_t>::max();
    int64_t end = 0;
    for (auto &v: variants)
    {
        if (v->getStart() < start) start = v->getStart();
        if (v->getEnd() > end) end = v->getEnd();
    }

    SimpleInterval combined_interval(contig, start, end);
    ReferenceContext combined_reference_context(cached_ref_seq_,
        combined_interval);
    combined_reference_context.setWindow(reference_window_padding,
        reference_window_padding);
    return combined_reference_context;
}
