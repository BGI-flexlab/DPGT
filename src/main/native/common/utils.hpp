#ifndef DPGT_UTILS_HPP
#define DPGT_UTILS_HPP


#include <string>
#include "htslib/kstring.h"

class Utils {
public:
    static void validateArg(bool condition, const std::string &msg);

    static
    std::string repeatJoin(const std::string &u, const std::string &s, int n);

    static void krepeatPuts(kstring_t *o, const std::string &u,
        const std::string &s, int n);
};


template <typename PtrT>
struct PtrLess {
    bool operator()(PtrT lhs, PtrT rhs) const
    {
        return *lhs < *rhs;
    }
};


template <typename T>
struct ArrayBuff {
    int l; // length
    int m; // capacity
    T *data;
};

#endif  // DPGT_UTILS_HPP
