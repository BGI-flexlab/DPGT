#include "math_utils.hpp"
#include "common/log10_factorial_cache.hpp"


JacobianLog10Table::JacobianLog10Table()
{
    int n_step = MAX_TOLERANCE * INV_STEP;
    cache = (double *)malloc(n_step*sizeof(double));
    for (int k = 0; k < n_step; ++k)
    {
        cache[k] = log10(1.0 + pow(10.0, -k*TABLE_STEP));
    }
}

JacobianLog10Table::~JacobianLog10Table()
{
    if (cache) free(cache);
}


Pow10Table::Pow10Table() {
    cache = (double *)aligned_alloc(8, POW10TABLE_SIZE*sizeof(double));
    for (int i = 0; i < POW10TABLE_SIZE; ++i) {
        cache[i] = pow(10.0f, -0.1*i);
    }
}

Pow10Table::~Pow10Table() {
    if (cache) free(cache);
}

double &Pow10Table::operator[](int index) const {
    return cache[index];
}


log10FactorialCache MathUtils::LOG10_FACTORIAL_CACHE = log10FactorialCache();

int MaxElementIndex(double *vec, long length, int start,
    int finish)
{
    if (vec == nullptr) return -1;  // fail
    if (start >= finish ) return -2;
    double max = -std::numeric_limits<double>::infinity();
    int max_i = 0;
    for (int i = start; i < finish; ++i)
    {
        if (vec[i] > max)
        {
            max = vec[i];
            max_i = i;
        }
    }
    return max_i;
}


int MinElement(int *vec, long length, int start, int finish)
{
    if (vec == nullptr) return -1;  // fail
    if (start >= finish ) return -2; // range error
    if (finish > length) return -3;  // finish out of range

    int min = std::numeric_limits<int>::max();
    for (int i = start; i < finish; ++i) {
        if (vec[i] < min) min = vec[i];
    }

    return min;
}


double Log10SumLog10(double x, double y)
{
    if (x < y)
    {
        std::swap(x, y);
    }

    return x + log10(1 + pow(10.0, y - x));
}

double Log10SumLog10(double x, double y, double z)
{
    if (x >= y && x >= z)  {  // x is max
        return x + log10(1 + pow(10.0, y - x) + pow(10.0, z - x));
    } else if (y >= z) {     // y is max
        return y + log10(1 + pow(10.0, x - y) + pow(10.0, z - y));
    } else {                 // z is max
        return z + log10(1 + pow(10.0, x - z) + pow(10.0, y - z));
    }
}

double Log10SumLog10(double *log10_values, long length)
{
    return Log10SumLog10(log10_values, length, 0, length);
}

double Log10SumLog10(double *log10_values, long length, long start)
{
    return Log10SumLog10(log10_values, length, start, length);
}

double Log10SumLog10(double *log10_values, long length, long start,
    long finish)
{
    if (start >= finish)
    {
        return -std::numeric_limits<double>::infinity();
    }

    int max_element_idx = MaxElementIndex(log10_values, length, start, finish);
    double max_value = log10_values[max_element_idx];

    if (max_value == -std::numeric_limits<double>::infinity())
    {
        return -std::numeric_limits<double>::infinity();
    }

    double sum = 1.0;
    for (int i = start; i < finish; ++i)
    {
        // double cur_val = log10_values[i];
        // if (i == max_element_idx ||
        //     cur_val == -std::numeric_limits<double>::infinity())
        // {
        //     continue;
        // }
        // double scaled_val = cur_val - max_value;
        // sum += pow(10.0f, scaled_val);
        sum += i == max_element_idx ? 0 : pow(10.0f, log10_values[i]-max_value);
    }

    // if (sum == std::numeric_limits<double>::quiet_NaN() ||
    //     sum == std::numeric_limits<double>::infinity())
    // {
    //     fprintf(stderr, "[%s] Error! log10 p: Values must be non-infinite and"
    //         " non-NAN", __func__);
    //     std::exit(1);
    // }

    return max_value + (sum != 1.0 ? log10(sum) : 0.0);
}


double ApproximateLog10SumLog10(double x, double y)
{
    static JacobianLog10Table jacobian_log10_table = JacobianLog10Table();
    if (x < y) std::swap(x, y);
    double diff = x - y;
    return x + (diff < jacobian_log10_table.MAX_TOLERANCE ?
        jacobian_log10_table.get(diff) : 0.0f);
}

double ApproximateLog10SumLog10(double x, double y, double z)
{
    return ApproximateLog10SumLog10(x, ApproximateLog10SumLog10(y, z));
}

double ApproximateLog10SumLog10(double *log10_values, long length)
{
    return ApproximateLog10SumLog10(log10_values, length, 0, length);
}

double ApproximateLog10SumLog10(double *log10_values, long length,
    long start)
{
    return ApproximateLog10SumLog10(log10_values, length, start, length);
}

double ApproximateLog10SumLog10(double *log10_values, long length,
    long start, long finish)
{
    static JacobianLog10Table jacobian_log10_table = JacobianLog10Table();
    if (start <= finish) return -std::numeric_limits<double>::infinity();
    int max_element_index = MaxElementIndex(log10_values, length, start, finish);
    double appro_sum = log10_values[max_element_index];

    for (int i = start; i < finish; ++i)
    {
        double cur_val = log10_values[i];
        if (i == max_element_index ||
            cur_val == -std::numeric_limits<double>::infinity())
        {
            continue;
        }
        double diff = appro_sum - cur_val;
        if (diff < jacobian_log10_table.MAX_TOLERANCE)
        {
            appro_sum += jacobian_log10_table.get(diff);
        }
    }
    return appro_sum;
}


double PLsumLog10(int a, int b)
{
    static const Pow10Table pow10_table = Pow10Table();
    if (a > b) {
        std::swap(a, b);
    }
    int s = b - a;
    if (s < POW10TABLE_SIZE) {
        return -0.1*a + log10(pow10_table[s]);
    } else {
        return -0.1*a + log10(pow(10.0f, -0.1*s));
    }
}

double PLsumLog10(int a, int b, int c)
{
    int pls[3] = {a, b, c};
    return PLsumLog10(pls, 3, 0, 3);
}


double PLsumLog10(int *pls, long length)
{
    return PLsumLog10(pls, length, 0, length);
}

double PLsumLog10(int *pls, long length, long start)
{
    return PLsumLog10(pls, length, start, length);
}

double PLsumLog10(int *pls, long length, long start, long finish)
{
    static const Pow10Table pow10_table = Pow10Table();
    if (start >= finish)
    {
        return -std::numeric_limits<double>::infinity();
    }

    int min_value = MinElement(pls, length, start, finish);

    double sum = 0.0;
    for (int i = start; i < finish; ++i)
    {
        int scaled_val = pls[i] - min_value;
        if (scaled_val < POW10TABLE_SIZE) {
            sum += pow10_table[scaled_val];
        } else {
            sum += pow(10.0f, -0.1*scaled_val);
        }

        // sum += pow(10.0f, log10_values[i]-max_value);
    }

    if (sum == std::numeric_limits<double>::quiet_NaN() ||
        sum == std::numeric_limits<double>::infinity())
    {
        fprintf(stderr, "[%s] Error! log10 p: Values must be non-infinite and"
            " non-NAN", __func__);
        std::exit(1);
    }

    // return max_value + (sum != 1.0 ? log10(sum) : 0.0);
    return -0.1*min_value + log10(sum);
}

