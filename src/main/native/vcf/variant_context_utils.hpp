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
#ifndef DPGT_VARIANT_CONTEX_UTILS_HPP
#define DPGT_VARIANT_CONTEX_UTILS_HPP


#include "common/simple_interval.hpp"
#include "vcf/allele.hpp"
#include "vcf/variant_context.hpp"
#include <vector>

/**
 * @brief 
 * 
 */
class VariantContextUtils {
public:

    /**
     * @brief Returns true if the context represents a MNP. If the context
     * supplied contains the GVCF NON_REF symbolic, allele, then the
     * determination is performed after discarding the NON_REF symbolic allele.
     * 
     * @param vc a Variant context
     * @return true if the context represents an unmixed MNP (i.e. all alleles
     * are non-symbolic and of the same length > 1), false otherwise
     */
    static bool isUnmixedMnpIgnoringNonRef(std::shared_ptr<VariantContext> &vc);

    /**
     * @brief get variant context range
     * 
     * @param vc a Variant context
     * @return SimpleInterval 
     */
    static SimpleInterval variantContextToInterval(
        std::shared_ptr<VariantContext> &vc);


    /**
     * @brief determine reference allele
     * 
     * @param vcs variants
     * @param loc position
     * @return Allele new reference allele
     */
    static Allele determineRefAllele(
        std::vector<std::shared_ptr<VariantContext>> &vcs,
        const SimpleInterval &loc);
    
    static const Allele &determineRefAllele(
        const Allele &ref1, const Allele &ref2);
    
    static bool variantMatchesLoc(
        std::shared_ptr<VariantContext> &vc, const SimpleInterval &loc);
};



#endif  // DPGT_VARIANT_CONTEX_UTILS_HPP
