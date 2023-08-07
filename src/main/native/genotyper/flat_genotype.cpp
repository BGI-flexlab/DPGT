#include "genotyper/flat_genotype.hpp"
#include "vcf/vcf_attribute.hpp"




void FlatGenotype::getString(const std::set<int> &format_key_indices,
    int ploidy, kstring_t *s) const
{   
    if (GT_ != nullptr) {
        GT_->getValueStr(s);
    } else {
        Utils::krepeatPuts(s, VCFConstants::MISSING_VALUE_v4.c_str(),
            VCFConstants::UNPHASED.c_str(), ploidy);
    }

    // Fixed the bug when there were no other fields except GT in FORMAT.
    // line:26 --last error
    if (format_key_indices.empty()) {
        return;
    }

    kputc(':', s);
    int n_miss = 0;
    VcfAttributeBase *a = nullptr;
    auto last = format_key_indices.end();
    --last;
    for (auto itr1 = format_key_indices.begin(); itr1 != last;
        ++itr1)
    {
        a = attributes_[*itr1];
        if (a) {
            if (n_miss > 0) {
                Utils::krepeatPuts(s,
                    VCFConstants::MISSING_VALUE_v4.c_str(),
                    ":", n_miss);
                n_miss = 0;
                kputc(':', s);
            }
            a->getValueStr(s);
            kputc(':', s);
        } else {
            ++n_miss;
        }
    }
    
    a = attributes_[*last];
    if (a) {
        if (n_miss > 0) {
            Utils::krepeatPuts(s,
                VCFConstants::MISSING_VALUE_v4.c_str(),
                ":", n_miss);
            n_miss = 0;
            kputc(':', s);
        }
        a->getValueStr(s);
    } else {
        s->s[--s->l] = 0;  // remove tail ':'
    }
}



