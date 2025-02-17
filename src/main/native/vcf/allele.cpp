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
#include <cctype>
#include <cstddef>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <utility>

#include "allele.hpp"
#include "nucleotide.hpp"


const std::string Allele::EMPTY_ALLELE_BASES = "";

const std::string Allele::NO_CALL_STRING = ".";

const std::string Allele::SPAN_DEL_STRING = "*";

const std::string Allele::NON_REF_STRING = "<NON_REF>";

const std::string Allele::UNSPECIFIED_ALTERNATE_ALLELE_STRING = "<*>";


Allele::Allele(const std::string &bases_, bool is_ref_) {
    is_null = false;
    if (wouldBeNullAllele(bases_)) {
        throw std::invalid_argument("Null alleles are not supported");
    }

    if (wouldBeNoCallAllele(bases_)) {
        this->bases = EMPTY_ALLELE_BASES;
        this->is_no_call = true;
        if (is_ref_) throw std::invalid_argument(
            "Cannot tag a NoCall allele as the reference allele");
        return;
    }

    if (wouldBeSymbolicAllele(bases_)) {
        this->is_symbolic = true;
        this->bases = bases_;
        std::transform(this->bases.begin(), this->bases.end(),
            this->bases.begin(), toupper);
        if (is_ref_) throw std::invalid_argument(
            "Cannot tag a symbolic allele as the reference allele");
        return;
    }

    this->is_ref = is_ref_;
    this->bases = bases_;
    std::transform(this->bases.begin(), this->bases.end(),
        this->bases.begin(), toupper);
    
    if (!acceptableAlleleBases(bases_, is_ref_))
        throw std::invalid_argument(
            "Unexcepted base in allele bases '" + bases_ + "'");
}

Allele::Allele(const Allele &allele, bool ignoreRefState)
{
    this->bases = allele.bases;
    this->is_null = allele.is_null;
    this->is_ref = ignoreRefState ? false : allele.is_ref;
    this->is_no_call = allele.is_no_call;
    this->is_symbolic = allele.is_symbolic;
}

Allele::Allele(const Allele &other) {
    this->bases = other.bases;
    this->is_null = other.is_null;
    this->is_ref = other.is_ref;
    this->is_no_call = other.is_no_call;
    this->is_symbolic = other.is_symbolic;
}

Allele::Allele(Allele &&other) noexcept {
    this->bases = std::move(other.bases);
    this->is_null = other.is_null;
    this->is_ref = other.is_ref;
    this->is_no_call = other.is_no_call;
    this->is_symbolic = other.is_symbolic;
}

Allele &Allele::operator=(const Allele &other) {
    if (&other != this) {
        this->bases = other.bases;
        this->is_null = other.is_null;
        this->is_ref = other.is_ref;
        this->is_no_call = other.is_no_call;
        this->is_symbolic = other.is_symbolic;
    }
    return *this;
}

Allele &Allele::operator=(Allele &&other) noexcept {
    if (&other != this) {
        this->bases = std::move(other.bases);
        this->is_null = other.is_null;
        this->is_ref = other.is_ref;
        this->is_no_call = other.is_no_call;
        this->is_symbolic = other.is_symbolic;
    }
    return *this;
}

// static const alleles used for optimization
const Allele Allele::NULL_ALLELE = Allele();
const Allele Allele::REF_A = Allele("A", true);
const Allele Allele::ALT_A = Allele("A", false);
const Allele Allele::REF_C = Allele("C", true);
const Allele Allele::ALT_C = Allele("C", false);
const Allele Allele::REF_G = Allele("G", true);
const Allele Allele::ALT_G = Allele("G", false);
const Allele Allele::REF_T = Allele("T", true);
const Allele Allele::ALT_T = Allele("T", false);
const Allele Allele::REF_N = Allele("N", true);
const Allele Allele::ALT_N = Allele("N", false);
const Allele Allele::SPAN_DEL = Allele(Allele::SPAN_DEL_STRING, false);
const Allele Allele::NO_CALL = Allele(Allele::NO_CALL_STRING, false);
const Allele Allele::NON_REF_ALLELE = Allele(Allele::NON_REF_STRING, false);
const Allele Allele::UNSPECIFIED_ALTERNATE_ALLELE = Allele(
    Allele::UNSPECIFIED_ALTERNATE_ALLELE_STRING, false);

static Allele OneRefBaseToAllele(char n)
{
    // "=ACMGRSVTWYHKDBN";
    static const Allele ref_allele_array[] ={
        Allele::NULL_ALLELE,
        Allele::REF_A,
        Allele::REF_C,
        Allele::NULL_ALLELE,
        Allele::REF_G,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::REF_T,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::REF_N};
    return ref_allele_array[seq_nt16_table[(int)n]];
}

static Allele OneAltBaseToAllele(char n)
{
    // "=ACMGRSVTWYHKDBN";
    static const Allele alt_allele_array[] ={
        Allele::NULL_ALLELE,
        Allele::ALT_A,
        Allele::ALT_C,
        Allele::NULL_ALLELE,
        Allele::ALT_G,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::ALT_T,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::NULL_ALLELE,
        Allele::ALT_N};
    return alt_allele_array[seq_nt16_table[(int)n]];
}


Allele Allele::create(const std::string &bases_, bool is_ref_) {
    if (bases_.empty()) {
        throw std::invalid_argument("bases can not be empty! use Allele() to "
            "create a empty allele!");
    }

    if (bases_.size() == 1) {
        // optimization for bases of size == 1.
        // optimize htsjdk java code by use lookup table to find allele of base
        switch (bases_[0]) {
            case '.':
                if (is_ref_) throw std::invalid_argument("Cannot tag a NoCall "
                    "allele as the reference allele");
                return NO_CALL;
            case '*':
                if (is_ref_) throw std::invalid_argument("Cannot tag a spanning"
                    " deletion allele as the reference allele");
                return SPAN_DEL;
            default:
                Allele allele;
                if (is_ref_) {
                    allele = OneRefBaseToAllele(bases_[0]);
                } else {
                    allele = OneAltBaseToAllele(bases_[0]);
                }
                if (allele.bases.empty()) {
                    throw std::invalid_argument("Illegal bases "
                        "'" + bases_ + "'");
                }
                return allele;
        }
    } else {
        return Allele(bases_, is_ref_);
    }
}

Allele Allele::create(char base, bool is_ref) {
    char bases_[1] = {base};
    return create(bases_, is_ref);
}

Allele Allele::create(const Allele &allele, bool ignoreRefState) {
    return Allele(allele, ignoreRefState);
}


Allele Allele::extend(const Allele &left, const std::string &rigth) {
    if (left.is_symbolic) throw std::invalid_argument("Cannot extend a "
        "symbolic allele");
    return create(left.bases+rigth, left.is_ref);
}


bool Allele::wouldBeNullAllele(const std::string &bases_) {
    return (bases_.size() == 1 && bases_[0] == VCFConstants::NULL_ALLELE)
        || bases_.size() == 0;
}

bool Allele::wouldBeStarAllele(const std::string &bases_) {
    return bases_.size() == 1 &&
        bases_[0] == VCFConstants::SPANNING_DELETION_ALLELE;
}

bool Allele::wouldBeNoCallAllele(const std::string &bases_) {
    return bases_.size() == 1 && bases_[0] == VCFConstants::NO_CALL_ALLELE;
}

bool Allele::wouldBeSymbolicAllele(const std::string &bases_) {
    if (bases_.size() <= 1) return false;
    return bases_.front() == SYMBOLIC_ALLELE_START ||
        bases_.back() == SYMBOLIC_ALLELE_END ||
        wouldBeBreakpoint(bases_) ||
        wouldBeSingleBreakend(bases_);
}

bool Allele::wouldBeBreakpoint(const std::string &bases_) {
    if (bases_.size() <= 1) {
        return false;
    }
    for (size_t i = 0; i < bases_.size(); i++) {
        if (bases_[i] == BREAKEND_EXTENDING_LEFT ||
            bases_[i] == BREAKEND_EXTENDING_RIGHT)
        {
            return true;
        }
    }
    return false;
}

bool Allele::wouldBeSingleBreakend(const std::string &bases_) {
    if (bases_.size() <= 1) return false;
    return bases_.front() == SINGLE_BREAKEND_INDICATOR ||
        bases_.back() == SINGLE_BREAKEND_INDICATOR;
}

bool Allele::acceptableAlleleBasesSimp(
    const std::string &bases_, bool is_ref_)
{
    if (wouldBeStarAllele(bases_)) return !is_ref_;
    for (auto const &b: bases_) {
        if(!nucleotide::IsStandardNt(b)) return false;
    }
    return true;
}

bool Allele::acceptableAlleleBases(
    const std::string &bases_, bool is_ref_)
{
    if (wouldBeNullAllele(bases_)) return false;

    if (wouldBeNoCallAllele(bases_) ||
        wouldBeSymbolicAllele(bases_)) return true;
    
    if (wouldBeStarAllele(bases_)) return !is_ref_;

    for (auto const &b: bases_) {
        // optimize htsjdk java code by use lookup table check if 
        // base is in "ACGTNacgtn"
        if(!nucleotide::IsStandardNtIncludeNn(b)) return false;
    }

    return true;
}

bool Allele::acceptableAlleleBases(const std::string &bases_) {
    return acceptableAlleleBases(bases_, true);
}

template<typename Iter>
Allele Allele::getMatchingAllele(Iter it, Iter end, const std::string &test) {
    for (; it != end; ++it) {
        if (it->basesMatch(test)) return *it;
    }

    if (wouldBeNoCallAllele(test))
        return NO_CALL;
    else
        return Allele();
}


int Allele::compareTo(const Allele &other) const {
    if (is_ref && !other.is_ref) {
        return -1;
    } else if (!is_ref && other.is_ref) {
        return 1;
    } else if (isNonRefAllele() && !other.isNonRefAllele()) {
        return 1;
    } else if (!isNonRefAllele() && other.isNonRefAllele()) {
        return -1;
    } else if (is_symbolic && !other.is_symbolic) {
        return 1;
    } else if (!is_symbolic && other.is_symbolic) {
        return -1;
    } else if (is_no_call && !other.is_no_call) {
        return 1;
    } else if (!is_no_call && other.is_no_call) {
        return -1;
    } else {
        return strcmp(getBasesString().c_str(), other.getBasesString().c_str());
    }
}


bool Allele::oneIsPrefixOfOther(const Allele &a1, const Allele &a2) {
    if (a2.length() >= a1.length()) {
        return firstIsPrefixOfSecond(a1, a2);
    } else {
        return firstIsPrefixOfSecond(a2, a1);
    }
}


bool Allele::firstIsPrefixOfSecond(const Allele &a1, const Allele &a2) {
    std::string a1_string = a1.getBasesString();
    return a2.getBasesString().substr(0, a1_string.size()) == a1_string;
}

