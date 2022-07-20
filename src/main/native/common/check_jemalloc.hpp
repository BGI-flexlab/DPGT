#ifndef DPGT_CHECK_JEMALLOC_HPP
#define DPGT_CHECK_JEMALLOC_HPP

#include <stddef.h>

extern "C"
{
    // weak symbol: resolved at runtime by the linker if we are using jemalloc, nullptr otherwise
    int mallctl(const char *name, void *oldp, size_t *oldlenp, void *newp, size_t newlen) __attribute__((weak));
}


bool check_jemalloc();


#endif  // DPGT_CHECK_JEMALLOC_HPP
