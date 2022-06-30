#ifndef DPGT_INT_TO_DOUBLE_CACHE_HPP
#define DPGT_INT_TO_DOUBLE_CACHE_HPP

#include <iostream>
#include <cstdlib>


/**
 * @brief abstract base class for int to double cache
 */
class IntToDoubleCache {
private:
    double *cache_ = nullptr;
    int size_ = 0;

public:
    IntToDoubleCache() = default;

    IntToDoubleCache(const IntToDoubleCache &other);
    IntToDoubleCache(IntToDoubleCache &&other) noexcept;

    IntToDoubleCache &operator=(const IntToDoubleCache &other);
    IntToDoubleCache &operator=(IntToDoubleCache &&other) noexcept;

    virtual ~IntToDoubleCache() {
        if (cache_) free(cache_);
    }
    
    virtual int maxSize() = 0;
    virtual double compute(int n) = 0;

    double get(int i);  // not thread safe

    void expandCache(int new_capacity);  // not thread safe

    int size() const {
        return size_;
    }
};



#endif  // DPGT_LOG10_FACTORIAL_CACHE_HPP
