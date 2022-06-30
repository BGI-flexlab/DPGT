#ifndef GLAT_SIMPLE_MATRIX_HPP
#define GLAT_SIMPLE_MATRIX_HPP


#include "utils.hpp"
#include <memory>
#include <cstdint>
#include <iostream>
#include <vector>


/**
 * A simple 2D matrix.
 * Datas are stored in 'row major'
 */
template<typename T>
class SimpleMatrix
{
public:
    enum class dimension: uint8_t {ROW, COLUMN};   

    SimpleMatrix() = default;
    
    SimpleMatrix(int64_t num_rows, int64_t num_cols):
        n_rows_(num_rows), n_cols_(num_cols)
    {
        if (num_rows < 0 || num_cols < 0)
        {
            std::cerr << "[SimpleMatrix] Error! Number of rows and cols must "
                << ">= 0, but the input is: " << num_rows << ", " << num_cols
                << std::endl;
            std::exit(1);
        }
        size_ = num_rows * num_cols;
        if (size_ > 0)
        {
            data_ = new T[size_];
        }
    }

    SimpleMatrix(T *input_data, int64_t num_rows, int64_t num_cols):
        n_rows_(num_rows), n_cols_(num_cols)
    {
        if (num_rows < 0 || num_cols < 0)
        {
            std::cerr << "[SimpleMatrix] Error! Number of rows and cols must "
                << ">= 0, but the input is: " << num_rows << ", " << num_cols
                << std::endl;
            std::exit(1);
        }
        size_ = num_rows * num_cols;
        if (size_ > 0)
        {
            data_ = new T[size_];
        }
        for (int64_t i = 0; i < size_; ++i)
        {
            data_[i] = input_data[i];
        }
    }

    SimpleMatrix(const SimpleMatrix &other)
    {
        n_rows_ = other.n_rows_;
        n_cols_ = other.n_cols_;
        size_ = other.size_;
        if (other.data_)
        {
            data_ = new T[size_];
            for (int64_t i = 0; i < size_; ++i)
            {
                data_[i] = other.data_[i];
            }
        } else
        {
            data_ = nullptr;
        }
    }

    SimpleMatrix(SimpleMatrix &&other) noexcept
    {
        n_rows_ = other.n_rows_;
        n_cols_ = other.n_cols_;
        size_ = other.size_;
        data_ = other.data_;
        other.data_ = nullptr;
    }

    SimpleMatrix &operator=(const SimpleMatrix &other)
    {
        if (this != &other)
        {
            n_rows_ = other.n_rows_;
            n_cols_ = other.n_cols_;
            size_ = other.size_;
            if (data_) delete [] data_;
            if (other.data_)
            {
                data_ = new T[size_];
                for (int64_t i = 0; i < size_; ++i)
                {
                    data_[i] = other.data_[i];
                }
            } else
            {
                data_ = nullptr;
            }
        }
        return *this;
    }

    SimpleMatrix &operator=(SimpleMatrix &&other) noexcept
    {
        if (this != &other)
        {
            n_rows_ = other.n_rows_;
            n_cols_ = other.n_cols_;
            size_ = other.size_;
            if (data_) delete [] data_;
            data_ = other.data_;
            other.data_ = nullptr;
        }
        return *this;
    }

    ~SimpleMatrix()
    {
        if (data_)
        {
            delete [] data_;
        }
    }


    /**
     * Access matrix value by index.
     */
    T &operator()(int64_t row_index, int64_t col_index)
    {
        row_index = row_index > 0 ? row_index : -row_index;
        col_index = col_index > 0 ? col_index : -col_index;
        if (row_index >= n_rows_ || col_index >= n_cols_)
        {
            std::cerr << "[SimpleMatrix] Error! Row index or column index "
                "out of range!" << std::endl;
            std::exit(1);
        }
        return data_[row_index * n_cols_ + col_index];
    }

    /**
     * Access matrix value by index.
     */
    const T &operator()(int64_t row_index, int64_t col_index) const
    {
        row_index = row_index > 0 ? row_index : -row_index;
        col_index = col_index > 0 ? col_index : -col_index;
        if (row_index >= n_rows_ || col_index >= n_cols_)
        {
            std::cerr << "[SimpleMatrix] Error! Row index or column index "
                "out of range!" << std::endl;
            std::exit(1);
        }
        return data_[row_index * n_cols_ + col_index];
    }

    /**
     * @brief get row data by row index
     */
    T *row(int64_t row_index) {
        row_index = row_index > 0 ? row_index : -row_index;
        if (row_index >= n_rows_)
        {
            std::cerr << "[SimpleMatrix] Error! Row index out of range!"
                << std::endl;
            std::exit(1);
        }
        return data_ + row_index * n_cols_;
    }

    /**
     * Returns a one dimension vector of this object's data in column major order
     */
    std::vector<T> ColumnMajorData() const
    {
        std::vector<T> result;
        for (int64_t j = 0; j < n_cols_; ++j)
        {
            for (int64_t i = 0; i < n_rows_; ++i)
            {
                result.push_back((*this)(i, j));
            }
        }
        return result;
    }

    /**
     * @brief fill this matrix with input value
     * 
     * @param val input value
     */
    void fill(const T &val) {
        for (int i = 0; i < size_; ++i)
            data_[i] = val;
    }

    /**
     * @brief fill subset of this matrix with input value
     * 
     * @param row_start inclusive start index of row
     * @param row_end exclusive end index of row
     * @param col_start inclusive start index of column
     * @param col_end exclusive end index of column
     * @param val input value
     */
    void fill(int64_t row_start, int64_t row_end,
        int64_t col_start, int64_t col_end, const T &val)
    {
        Utils::validateArg(row_start >= 0 && row_start < n_rows_,
            "[SimpleMatrix] Error! row_start out of range");
        Utils::validateArg(row_end > 0 && row_end <= n_rows_,
            "[SimpleMatrix] Error! row_end out of range");
        Utils::validateArg(col_start >= 0 && col_start < n_cols_,
            "[SimpleMatrix] Error! row_end out of range");
        Utils::validateArg(col_end > 0 && col_end <= n_cols_,
            "[SimpleMatrix] Error! col_end out of range");
        for (int i = row_start; i < row_end; ++i) {
            const int first_row_index = i * n_cols_;
            for (int j = col_start; j < col_end; ++j) {
                data_[first_row_index+j] = val;
            }
        }
    }

    /**
     * Number of rows of the matrix.
     */
    const int64_t &rows() const
    {
        return n_rows_;
    }

    /**
     * Number of columns of the matrix
     */
    const int64_t &cols() const
    {
        return n_cols_;
    }

    const T *data() const
    {
        return data_;
    }

private:
    T *data_ = nullptr;
    int64_t n_rows_;  // number of rows
    int64_t n_cols_;  // number of cols
    int64_t size_;
};


template<typename T>
std::ostream &operator<<(std::ostream &os, const SimpleMatrix<T> &mat)
{
    for (int64_t i = 0; i < mat.rows(); ++i)
    {
        for (int64_t j = 0; j < mat.cols(); ++j)
        {
            os << mat(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}


#endif  // GLAT_SIMPLE_MATRIX_HPP

