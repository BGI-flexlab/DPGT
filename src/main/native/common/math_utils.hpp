#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP

#include <limits>
#include <math.h>
#include <vector>
#include <memory>
#include <cmath>
#include <map>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "log10_factorial_cache.hpp"

class MathUtils {
public:
    static log10FactorialCache LOG10_FACTORIAL_CACHE;

    static double log10Factorial(int n) {
        return LOG10_FACTORIAL_CACHE.get(n);
    }
};


// common eigen types that we will often use
template <typename T>
using EigenMatrixMap =
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;

template <typename T>
using EigenArrayMap =
    Eigen::Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>>;

template <typename T>
using EigenVectorMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>;

template <typename T>
using EigenVectorArrayMap = Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>;

template <typename T>
using ConstEigenMatrixMap =
    Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;

template <typename T>
using ConstEigenArrayMap =
    Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>>;

template <typename T>
using ConstEigenVectorMap =
    Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>;

template <typename T>
using ConstEigenVectorArrayMap =
    Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>>;

using EigenOuterStride = Eigen::OuterStride<Eigen::Dynamic>;
using EigenInnerStride = Eigen::InnerStride<Eigen::Dynamic>;
using EigenStride = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using EigenOuterStridedMatrixMap = Eigen::
    Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, 0, EigenOuterStride>;

template <typename T>
using EigenOuterStridedArrayMap = Eigen::
    Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>, 0, EigenOuterStride>;

template <typename T>
using ConstEigenOuterStridedMatrixMap = Eigen::Map<
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>,
    0, EigenOuterStride>;

template <typename T>
using ConstEigenOuterStridedArrayMap = Eigen::Map<
    const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>,
    0, EigenOuterStride>;

template <typename T>
using EigenStridedMatrixMap = Eigen::
    Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, 0, EigenStride>;

template <typename T>
using EigenStridedArrayMap =
    Eigen::Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>, 0, EigenStride>;

template <typename T>
using ConstEigenStridedMatrixMap = Eigen::
    Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, 0, EigenStride>;

template <typename T>
using ConstEigenStridedArrayMap = Eigen::
    Map<const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>, 0, EigenStride>;

// 1-d array of dynamic length
template <typename T>
using EigenArrayXt = Eigen::Array<T, Eigen::Dynamic, 1>;

using EigenArrayXf = Eigen::ArrayXf;
using EigenArrayXd = Eigen::ArrayXd;
using EigenArrayXi = Eigen::ArrayXi;
using EigenArrayXb = EigenArrayXt<bool>;


// 2-d array of dynamic size
template <typename T>
using EigenArrayXXt = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>;

using EigenArrayXXf = Eigen::ArrayXXf;
using EigenArrayXXd = Eigen::ArrayXXd;
using EigenArrayXXi = Eigen::ArrayXXi;
using EigenArrayXXb = EigenArrayXXt<bool>;


/**
 * Returns a column major one dimension data given a row major one dimension
 * data.
 * @param data row major one dimension array which represent a 2D matrix.
 * @param x number of rows.
 * @param y number of columns.
 * @return column major one dimension data.
 */
template <typename T>
std::vector<T> ColumnMajorMatrixData(const T *data, int64_t x, int64_t y)
{
    std::vector<T> result;
    result.reserve(x*y);
    for (int64_t j = 0; j < y; ++j)
    {
        for (int64_t i = 0; i < x; ++i)
        {
            result.push_back(data[i * y + j]);
        }
    }
    return result;
}


/**
 * Composes and array of doubles given the start, end (enclosed*) and the step.
 * 
 * For example Doubles(1, 10, 1.0) results in the sequence
 * [1.0, 2.0, 3.0, ..., 9.0, 10.0]
 * 
 * It also works with a negative step as long as end < start.
 * For example Doubles(3, 1, -.5)} results in [3.0, 2.5, 2.0, 1.5, 1.0]
 * 
 * The 'end' values might be included if there is a N >= 1 such that
 * N * step + start == end. There is a difference "tolerances" three
 * orders of magnitude below thet step absolut value's so if the step is 1.
 * and the last value is under 0.001 absolute difference with 'end',
 * the 'end' value is used instead.
 *
 * A sequence request where the difference between the start and end is less
 * than the 'tolerance' described above results in a single value sequence
 * equal to the start.
 * 
 * @param start the first values of the sequence.
 * @param limit the limit of the sequence
 * @param step the step increment (or decrement) between elements of the sequence.
 * @return sequence of doubles.
 */
std::vector<double> Doubles(double start, double limit, double step);


/**
 * phred score calculations
 */

const double kLog10OneHalf = -0.3010299956639812;  // log10(0.5);

class JacobianLog10Table {
public:
    JacobianLog10Table();
    ~JacobianLog10Table();
    // if log(a) - log(b) > MAX_TOLERANCE, b is effectively treated as zero in
    // approximateLogSumLog
    // MAX_TOLERANCE = 8.0 introduces an error of at most one part in 10^8 in sums
    double MAX_TOLERANCE = 8.0;

    double get(double diff) const
    {
        int index = round(diff * INV_STEP);
        return cache[index];
    }

private:
    double TABLE_STEP = 0.0001;
    double INV_STEP = 10000;  // 1/TABLE_STEP
    double *cache = nullptr;
};


#define POW10TABLE_SIZE 2048
#define POW10TABLE_SIZE_F 204.8f

class Pow10Table {
public:
    Pow10Table();
    ~Pow10Table();

    double &operator[](int index) const;

private:
    double *cache = nullptr;
};


int MaxElementIndex(double *vec, long length, int start,
    int finish);

int MinElement(int *vec, long length, int start, int finish);

double Log10SumLog10(double x, double y);

double Log10SumLog10(double x, double y, double z);

double Log10SumLog10(double *log10_values, long length);

double Log10SumLog10(double *log10_values, long length, long start);

double Log10SumLog10(double *log10_values, long length, long start,
    long finish);


double ApproximateLog10SumLog10(double x, double y);

double ApproximateLog10SumLog10(double x, double y, double z);

double ApproximateLog10SumLog10(double *log10_values, long length);

double ApproximateLog10SumLog10(double *log10_values, long length,
    long start);

double ApproximateLog10SumLog10(double *log10_values, long length,
    long start, long finish);


double PLsumLog10(int a, int b);
double PLsumLog10(int a, int b, int c);
double PLsumLog10(int *pls, long length);
double PLsumLog10(int *pls, long length, long start);
double PLsumLog10(int *pls, long length, long start, long finish);


inline
double Log10OneMinusPow10(double x)
{
    return log10(1 - pow(10, x));
}

#define LOG10_2 0.3010299956639812
#define LOG10_INFORMATIVE_THREHOLD 0.20411998265592482  // log10(1.6)
#define LOG10_3 0.47712125471966244
#define LOG10_6 0.7781512503836436


/**
 * Compare two values. input values must have operator< and operator== defined.
 * @param x the first value to compare.
 * @param y the second value to compare.
 * @return 0 if x ==y; -1 if x < y; 1 if x > y.
 */
template<typename T>
int Compare(const T &x, const T &y)
{
    return (x < y) ? -1 : ((x == y) ? 0 : 1);
}

template <typename T>
double RootMeanSquare(const EigenArrayXt<T> &v)
{
    return sqrt( (v*v).sum() / static_cast<double>(v.size()) );
}

#define MAX_LOG10_RATIO 255

inline
double Log10Ratio(int x, int y)
{
    if (x == 0 && y == 0) return 0.0;
    if (x < y) std::swap(x, y);
    if (y == 0) return MAX_LOG10_RATIO;
    double res = log10( static_cast<double>(x) / y );
    res = res > MAX_LOG10_RATIO ? MAX_LOG10_RATIO : res;
    return res;
}

template <typename T, typename T1>
T1 VecRound(const T &vec)
{
    T1 res(vec.size());
    for (int i = 0; i < static_cast<int>(vec.size()); ++i)
    {
        res[i] = static_cast<int>(std::round(vec[i]));
    }
    return res;
}

template <typename T>
double Sum(const std::vector<T> &vec)
{
    double sum = 0.0;
    for (auto const &i: vec)
    {
        sum += i;
    }
    return sum;
}

template <typename T>
double Mean(const std::vector<T> &vec)
{
    assert(!vec.empty());
    return Sum(vec) / vec.size();
}

template <typename T>
T Max(const std::vector<T> &vec)
{
    assert(!vec.empty());
    T max = -std::numeric_limits<T>::max();
    for (auto const &i: vec) {
        if (i > max) max = i;
    }
    return max;
}

template <typename T>
T Min(const std::vector<T> &vec)
{
    assert(!vec.empty());
    T min = std::numeric_limits<T>::max();
    for (auto const &i: vec) {
        if (i < min) min = i;
    }
    return min;
}


#endif  // MATH_UTILS_HPP
