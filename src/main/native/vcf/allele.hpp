#ifndef DPGT_ALLELE_HPP
#define DPGT_ALLELE_HPP

#include <algorithm>
#include <cstring>
#include <string>
#include "vcf_constants.hpp"


class Allele {
private:
    static const std::string EMPTY_ALLELE_BASES;
    static const char SINGLE_BREAKEND_INDICATOR = '.';
    static const char BREAKEND_EXTENDING_RIGHT = '[';
    static const char BREAKEND_EXTENDING_LEFT = ']';
    static const char SYMBOLIC_ALLELE_START = '<';
    static const char SYMBOLIC_ALLELE_END = '>';

    std::string bases = "";

    bool is_null = true;  // null allele constructed by default constructor

    bool is_ref = false;
    bool is_no_call = false;
    bool is_symbolic = false;

private:
    /**
     * Assume that input bases is not null allele, NO_CALL allele or
     * Symbolic allele. Use this function in constructor to avoid check
     * null allele, NO_CALL allele or Symbolic allele(would*) repeatly.
     * @param bases_ bases representing an allele
     * @param is_ref_ true if a reference allele
     * @return true if the bases represent well formated allele
     */
    static bool acceptableAlleleBasesSimp(
        const std::string &bases_, bool is_ref_);

public:
    static const std::string NO_CALL_STRING;
    
    static const std::string SPAN_DEL_STRING;

    static const std::string NON_REF_STRING;

    static const std::string UNSPECIFIED_ALTERNATE_ALLELE_STRING;

    Allele() = default;

    Allele(const std::string &bases_, bool is_ref_);

    /**
     * Creates a new allele based on the provided one. Ref state will be
     * copied unless ignoreRefState is true
     * (in which case the returned allele will be non-Ref).
     *
     * This method is efficient because it can skip the validation of the
     * bases (since the original allele was already validated).
     *
     * @param allele  the allele from which to copy the bases
     * @param ignoreRefState  should we ignore the reference state of the
     * input allele and use the default ref state?
     */
    Allele(const Allele &allele, bool ignoreRefState);

    Allele(const Allele &other);

    Allele(Allele &&other) noexcept;

    Allele &operator=(const Allele &other);

    Allele &operator=(Allele &&other) noexcept;

    ~Allele() {}

    bool operator<(const Allele &other) const {
        return compareTo(other) < 0;
    }

    bool operator>(const Allele &other) const {
        return compareTo(other) > 0;
    }

    bool operator==(const Allele &other) const {
        return compareTo(other) == 0;
    }

    // static const alleles used for optimization
    static const Allele NULL_ALLELE;
    static const Allele REF_A;
    static const Allele ALT_A;
    static const Allele REF_C;
    static const Allele ALT_C;
    static const Allele REF_G;
    static const Allele ALT_G;
    static const Allele REF_T;
    static const Allele ALT_T;
    static const Allele REF_N;
    static const Allele ALT_N;
    static const Allele SPAN_DEL;
    static const Allele NO_CALL;
    static const Allele NON_REF_ALLELE;
    static const Allele UNSPECIFIED_ALTERNATE_ALLELE;

    // omit symbolic Structure Vatriant defines

    /**
     * Create a new Allele that includes bases and if tagged as the reference
     * allele if isRef == true.  If bases == '-', a Null allele is created.
     * If bases ==  '.', a no call Allele is created.
     * If bases ==  '*', a spanning deletions Allele is created.
     *
     * @param bases the DNA sequence of this variation, '-', '.', or '*'
     * @param isRef should we make this a reference allele?
     */
    static Allele create(const std::string &bases_, bool is_ref_ = false);

    static Allele create(char base, bool is_ref_ = false);

    static Allele create(const Allele &allele, bool ignoreRefState);

    static Allele extend(const Allele &left, const std::string &rigth);

    /**
     * @brief
     * null allele: bases is "-" or bases is empty
     * @param bases_ bases representing an allele
     * @return true if the bases represent the null allele
     */
    static bool wouldBeNullAllele(const std::string &bases_);

    /**
     * @brief
     * SPAN_DEL allele: bases is "*".
     * @param bases_ bases representing an allele
     * @return true if the bases represent the SPAN_DEL allele
     */
    static bool wouldBeStarAllele(const std::string &bases_);

    /**
     * @brief 
     * NO_CALL allele: bases is "."
     * @param bases_ bases representing an allele
     * @return true if the bases represent the NO_CALL allele
     */
    static bool wouldBeNoCallAllele(const std::string &bases_);

    /**
     * @brief 
     * @param bases_ bases representing an allele
     * @return true if the bases represent a symbolic allele,
     * including breakpoints and breakends
     */
    static bool wouldBeSymbolicAllele(const std::string &bases_);

    /**
     * @brief 
     * @param bases_ bases representing an allele
     * @return true  if the bases represent a symbolic allele in breakpoint
     * notation, (ex: G]17:198982] or ]13:123456]T )
     */
    static bool wouldBeBreakpoint(const std::string &bases_);

    /**
     * @brief 
     * @param bases_ bases representing an allele
     * @return true if the bases represent a symbolic allele in single
     * breakend notation (ex: .A or A. )
     */
    static bool wouldBeSingleBreakend(const std::string &bases_);

    /**
     * @param bases_ bases representing an allele
     * @param is_ref_ true if a reference allele
     * @return true if the bases represent well formated allele
     */
    static bool acceptableAlleleBases(
        const std::string &bases_, bool is_ref_);
    
    /**
     * @param bases_ bases representing a reference allele
     * @return true if the bases represent well formated allele
     */
    static bool acceptableAlleleBases(const std::string &bases_);


    // accessor functions

    bool isNull() const {
        return is_null;
    }

    /**
     * @return true if this allele is a NO_CALL allele
     */
    bool isNoCall() const {
        return is_no_call;
    }

    /**
     * @return true if this allele is not a NO_CALL allele
     */
    bool isCalled() const {
        return !is_no_call;
    }

    /**
     * @return true if this allele is a reference allele
     */
    bool isReference() const {
        return is_ref;
    }

    /**
     * @return true if this allele is not a reference allele 
     *
     * Note that this is not the same as <NON_REF> or <*> allele.
     */
    bool isNonReference() const {
        return !is_ref;
    }

    /**
     * @return true if this allele is a symbolic allele(<.*>, breakpoints,
     * breakends)
     */
    bool isSymbolic() const {
        return is_symbolic;
    }

    /**
     * @return true if this allele is a breakpoint ( ex: G]17:198982]
     * or ]13:123456]T )
     */
    bool isBreakpoint() const {
        return wouldBeBreakpoint(bases);
    }

    /**
     * @return true if this Allele is a single breakend (ex: .A or A.)
     */
    bool isBreakend() const {
        return wouldBeSingleBreakend(bases);
    }

    /**
     * @return std::string a nice string representation of this object
     * append a "*" to reference allele.
     */
    std::string toString() const {
        return (is_no_call ? NO_CALL_STRING : bases) + (is_ref ? "*" : "");
    }

    /**
     * @return std::string empty string if this allele is symbolic allele or
     * NO_CALL allele; otherwise return bases
     */
    std::string getBases() const {
        return is_symbolic ? EMPTY_ALLELE_BASES : bases;
    }

    /**
     * @brief different from getBases, this function return "." if this 
     * allele is NO_CALL.
     * @return std::string 
     */
    std::string getBasesString() const {
        return is_no_call ? NO_CALL_STRING : getBases();
    }
    
    /**
     * @return std::string the bases.
     */
    std::string getDisplayString() const {
        return bases;
    }

    /**
     * @return true if this allele is identical to the other allele
     */
    bool equals(const Allele &other) const {
        return this->is_ref == other.is_ref &&
            this->is_no_call == other.is_no_call &&
            this->is_symbolic == other.is_symbolic &&
            this->bases == other.bases;
    }

    /**
     * @param test bases to test against
     * @return true if this allele contains the same bases as test, regardless 
     * of its reference status; handles Null and NO_CALL alleles.
     * Note that always return false if this allele is a symbolic allele.
     */
    bool basesMatch(std::string test) {
        std::transform(test.begin(), test.end(), test.begin(), toupper);
        return !is_symbolic && test == bases;
    }

    /**
     * @param test allele to test against
     * @return true true if this allele contains the same bases as test, regardless 
     * of its reference status; handles Null and NO_CALL alleles.
     * Note that always return false if this allele is a symbolic allele.
     */
    bool basesMatch(const Allele &test) {
        return !is_symbolic && test.getBases() == bases;
    }


    int length() const {
        return is_symbolic ? 0 : bases.length();
    }

    template<typename Iter>
    static Allele getMatchingAllele(Iter it, Iter end, const std::string &test);


    /**
     * @param other allele to compare with
     * @return int < 0 if this allele is less than other,
     *             > 0 if this allele is greater than other,
     *             = 0 if this allele is equal to other.
     * reference allele is less than non reference allele
     * Note that compareTo not use is_no_call and is_symbolic to do comparison
     */
    int compareTo(const Allele &other) const;

    static bool oneIsPrefixOfOther(const Allele &a1, const Allele &a2);

    static bool firstIsPrefixOfSecond(const Allele &a1, const Allele &a2);


    /**
     * @return true if this is a <NON_REF> or <*> allele
     */
    bool isNonRefAllele() const {
        return this->equals(NON_REF_ALLELE) ||
            this->equals(UNSPECIFIED_ALTERNATE_ALLELE);
    }

};


#endif  // DPGT_ALLELE_HPP
