#include "vcf/variant_context_utils.hpp"
#include "common/simple_interval.hpp"
#include "vcf/allele.hpp"
#include "vcf/variant_context.hpp"
#include <vector>



bool VariantContextUtils::isUnmixedMnpIgnoringNonRef(
    std::shared_ptr<VariantContext> &vc)
{
    const std::vector<Allele> &alleles = vc->getAlleles();
    int length = vc->getReference().length();
    if (length < 2) return false;
    for (auto &a: alleles) {
        if (a.isSymbolic() && !a.isNonRefAllele()) return false;
        else if (!a.isSymbolic() && a.length() != length) return false;
    }
    return true;
}


SimpleInterval VariantContextUtils::variantContextToInterval(
    std::shared_ptr<VariantContext> &vc)
{
    return SimpleInterval(vc->getContig(), vc->getStart(), vc->getEnd());
}


Allele VariantContextUtils::determineRefAllele(
    std::vector<std::shared_ptr<VariantContext>> &vcs,
    const SimpleInterval &loc)
{
    Allele ref;
    for (auto &vc: vcs) {
        // variant.start == loc.start
        if (variantMatchesLoc(vc, loc)) {
            Allele my_ref = vc->getReference();
            ref = determineRefAllele(ref, my_ref);
        }
    }

    return ref;
}

const Allele &VariantContextUtils::determineRefAllele(
    const Allele &ref1, const Allele &ref2)
{
    if (ref1.isNull() || ref1.length() < ref2.length()) {
        return ref2;
    } else if (ref2.isNull() || ref2.length() < ref1.length()) {
        return ref1;
    } else if (ref1.length() == ref2.length() && !ref1.equals(ref2)) {
        std::cerr << "[VariantContextUtils::determineRefAllele] Error! "
            << "the input reference alleles do not appear to represent the same"
            << " position, " << ref1.toString() << ", " << ref2.toString()
            << std::endl;
        std::exit(1);
    } else {
        return ref1;
    }
}

bool VariantContextUtils::variantMatchesLoc(
    std::shared_ptr<VariantContext> &vc, const SimpleInterval &loc)
{
    return loc.IsNull() || vc->getStart() == loc.start;
}