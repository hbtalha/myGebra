#include "Matrix.h"

#include <iostream>
#include <assert.h>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <boost/format.hpp>

using namespace std;

namespace L_Algebra
{

Matrix::Matrix(std::initializer_list<std::initializer_list<Fraction>> values )
{
    rows_num = values.size();
    cols_num = 0;

    for(auto &row: values)
    {
        if(row.size() != 0)
        {
            cols_num = row.size();
            break;
        }
    }

    assert(cols_num > 0);

    data.reserve(rows_num * cols_num);

    for(auto &row: values)
    {
        if(row.size() == 0)
            for(size_t i = 0; i < cols_num; ++i)
                data.push_back(0);
        else
        {
            assert(row.size() == cols_num);
            std::copy(row.begin(), row.end(), std::back_inserter(data));
        }
    }
}

bool Matrix::has_one_row_zero() const
{
    bool has;

    for(size_t i = 0; i < rows_num; ++i)
    {
        has = true;
        for(size_t j = 0; j < cols_num; ++j)
            if(at(i,j) != 0)
            {
                has = false;
                break;
            }

        if(has)
            return true;
    }

    return false;
}

ostream& operator<<(ostream& os, const Matrix& mx)
{
    size_t width = 1;
    for(const auto element : mx.data)
    {
        auto w = element.to_string().size();
        if(width < w)
            width = w;
    }

    string w = "%" + to_string(width + 4) + "d";

    for (size_t i = 0; i < mx.rows(); i++)
    {
        for (size_t j = 0; j < mx.cols(); j++)
            os << boost::format(w.c_str()) %  mx.at(i, j);

        os << '\n';
    }

    return os;
}

// to print the diagonal
std::ostream& operator<<(std::ostream& os,  const std::vector<Fraction>& v)
{
    for (auto e: v)
        os << e << " ";

    return os;
}

Matrix Matrix::operator+(const Matrix& mx) const
{
    assert(rows_num == mx.rows_num && cols_num == mx.cols_num);

    Matrix addition(rows_num, cols_num);

    transform(data.begin(), data.end(), mx.data.begin(), addition.data.begin(), plus{});

    return addition;
}

Matrix Matrix::operator-(const Matrix& mx) const
{
    assert(rows_num == mx.rows_num && cols_num == mx.cols_num);

    Matrix subtraction(rows_num, cols_num);

    transform(data.begin(), data.end(), mx.data.begin(), subtraction.data.begin(), plus{});

    return subtraction;
}

Matrix Matrix::operator*(const Matrix& mx) const
{
    assert(cols_num == mx.rows_num);

    Matrix multiplication(rows_num, mx.cols_num);

    for(size_t i = 0; i < rows_num; ++i)
        for (size_t j = 0; j < mx.cols_num; ++j)
            for(size_t x = 0; x < cols_num; ++x)
                multiplication.at(i,j) += at(i, x) * mx.at(x, j);

    return multiplication;
}

Matrix& Matrix::operator*=(const Matrix& mx)
{
    assert(cols_num == mx.rows_num);

    return *this = (*this * mx);
}

Matrix& Matrix::operator-=(const Matrix& mx)
{
    assert(rows_num == mx.rows_num && cols_num == mx.cols_num);

    transform(data.begin(), data.end(), mx.data.begin(), data.begin(), minus{});

    return *this;
}

Matrix& Matrix::operator+=(const Matrix& mx)
{
    assert(rows_num == mx.rows_num && cols_num == mx.cols_num);

    transform(data.begin(), data.end(), mx.data.begin(), data.begin(), plus{});

    return *this;
}

Matrix operator*(const Matrix& mx, Fraction n)
{
    Matrix multiplication(mx.rows_num, mx.cols_num);

    for(size_t i = 0; i < mx.rows_num; ++i)
        for(size_t j = 0; j < mx.cols_num; ++j)
            multiplication.at(i, j) = mx.at(i, j) * n;

    return multiplication;
}

Matrix operator*(Fraction n, const Matrix& mx)
{
    Matrix multiplication(mx.rows_num, mx.cols_num);

    for(size_t i = 0; i < mx.rows_num; ++i)
        for(size_t j = 0; j < mx.cols_num; ++j)
            multiplication.at(i, j) = mx.at(i, j) * n;

    return multiplication;
}

Matrix& Matrix::operator*=(const Fraction& n)
{
    return *this = *this * n;
}

Matrix Matrix::operator/(const Fraction& n) const
{
    assert(n != 0);

    Matrix division(rows_num, cols_num);

    for(size_t i = 0; i < rows_num; ++i)
        for(size_t j = 0; j < cols_num; ++j)
            division.at(i, j) = at(i, j) / n;

    return division;
}

Matrix Matrix::Identity(int n)
{
    assert(n > 0);

    Matrix mx(n,n);

    for(int i = 0; i < n; ++i)
        mx.at(i, i) = 1;

    return mx;
}

Matrix Matrix::Constant(int r, int c, long long n)
{
    Matrix mx(r,c, n);

    return mx;
}

bool Matrix::is_identity() const
{
    if(! is_square())
        return false;

    return *this == Identity(cols_num);
}

bool Matrix::is_symmetric() const
{
    if(! is_square())
        return false;

    for(size_t i = 0; i < rows_num; ++i)
        for(size_t j = 0; j < cols_num; ++j)
            if(at(i,j) != at(j,i))
                return false;

    return true;
}

bool Matrix::is_skewSymmetric() const
{
    if(! is_square())
        return false;

    for(size_t i = 0; i < rows_num; ++i)
        for(size_t j = i+1; j < cols_num; ++j)
            if(at(i,j) != -at(j,i))
                return false;

    return true;
}

bool Matrix::is_diagonal() const
{
    if(! is_square())
        return false;

    for(size_t i = 0; i < rows_num; ++i)
        for(size_t j = 0; j < cols_num; ++j)
            if(i != j)
                if( at(i, j) != 0 )
                    return false;

    return true;
}

bool Matrix::is_zero() const
{
    return all_of( data.begin(), data.end(), [ ] (const auto& x)
    {
        return x == 0;
    } );
}

bool Matrix::is_constant() const
{
    return adjacent_find( data.begin(), data.end(), not_equal_to{} ) == data.end();
}

bool Matrix::is_orthogonal() const
{
    if(! is_square())
        return false;

    return (*this * transpose() == Identity(cols_num));
}

bool Matrix::is_invertible() const
{
    return this->determinant() != 0;
}

bool Matrix::is_linearly_dependent() const
{
    return this->determinant() == 0;
}

bool Matrix::is_linearly_independent() const
{
    return ! this->is_linearly_dependent();
}

bool Matrix::is_lowerTriangular() const
{
    if(! is_square())
        return false;

    for(size_t i = 0; i < rows_num; ++i)
        for(size_t j = i + 1; j < cols_num; ++j)
            if( at(i,j) )
                return false;

    return true;
}

bool Matrix::is_upperTriangular() const
{
    if(! is_square())
        return false;

    for(size_t i = 0; i < rows_num; ++i)
        for(size_t j = 0; j < i; ++j)
            if( at(i,j) )
                return false;

    return true;
}

bool Matrix::is_consistent( ) const
{
    Matrix mx1 = gaussJordanElimination();

    bool square = is_square();

    int num_non_zero_numbers = 0;
    for(size_t i = 0; i < rows_num; ++i)
    {
        if (square)
            for(size_t j = 0; j < cols_num; ++j)
            {
                if(mx1(i, j) != 0)
                    ++num_non_zero_numbers;
            }
        else
            for(size_t j = 0; j < cols_num - 1; ++j)
            {
                if(mx1(i, j) != 0)
                    ++num_non_zero_numbers;
            }

        if( ! square && num_non_zero_numbers == 0 && mx1(i, cols_num - 1) != 0)
            return false;

        if(num_non_zero_numbers > 1)
            return false;

        num_non_zero_numbers = 0;
    }

    return true;
}

Matrix Matrix::transpose() const
{
    Matrix trans(cols_num, rows_num);

    for(size_t i = 0; i < rows_num; ++i)
        for(size_t j = 0; j < cols_num; ++j)
            trans.at(j, i) = at(i, j);

    return trans;
}

Fraction Matrix::trace() const
{
    assert(is_square());

    Fraction tr;
    for(size_t i = 0; i < rows_num; ++i)
        tr += at(i,i);

    return tr;
}

size_t Matrix::rank() const
{
    Matrix mx = this->gaussJordanElimination();

    int rank = 0;

    for(size_t i = 0; i < rows_num; ++i)
        for(size_t j = 0; j < cols_num; ++j)
            if(mx(i, j) != 0)
            {
                ++rank;

                break;
            }

    return rank;
}

Fraction Matrix::determinant() const
{
    assert(is_square());

    if(is_zero())
        return {0};

    if(has_one_row_zero())
        return {0};

    if(rows_num == 1)
        return at(0,0);

    if(is_identity())
        return {1};

    if(is_constant())
        return {0};

    if(cols_num == 2)
        return at(0,0) * at(1,1) - at(0,1) * at(1,0);

    bool alternative_pivot_1_found;

    bool pivot_not_zero_found;

    bool number_not_zero_found;

    size_t row_with_alternative_pivot;

    size_t row_with_pivot_not_zero;

    size_t pivot_row = 0;
    size_t pivot_col = 0;

    Matrix mx(*this);
    vector<Fraction> row_mults;
    int sign = 1;

    while (pivot_row < (rows_num - 1))
    {
        alternative_pivot_1_found = mx.pivotEqualTo_one_Found ( pivot_row, pivot_col, row_with_alternative_pivot);

        pivot_not_zero_found = mx.pivotNot_zero_Found(pivot_row, pivot_col, row_with_pivot_not_zero);

        if (mx.at(pivot_row, pivot_col) != 1 && alternative_pivot_1_found )
        {
            mx.swapRows(pivot_row, row_with_alternative_pivot);

            sign *= (-1);
        }
        else if (mx.at(pivot_row, pivot_col) == 0 && pivot_not_zero_found )
        {
            mx.swapRows(pivot_row, row_with_pivot_not_zero);

            sign *= (-1);
        }

        size_t col_dif_zero;

        number_not_zero_found = mx.firstNumberNot_zero(pivot_row, col_dif_zero);

        if(number_not_zero_found)
        {
            if (mx.at(pivot_row, col_dif_zero) != 1)
            {
                row_mults.push_back(mx.at(pivot_row, col_dif_zero));

                mx.changePivotTo_one(pivot_row, mx.at(pivot_row, col_dif_zero));
            }
        }

        for (size_t i = pivot_row + 1; i < rows_num; ++i)
            mx.zeroOutTheColumn(i, pivot_row, mx.at(i, col_dif_zero));

        ++pivot_row;
        ++pivot_col;
    }

    Fraction det(sign);

    for(size_t i = 0; i < rows_num; ++i)
        det  *= mx.at(i,i);

    return accumulate(row_mults.begin(), row_mults.end(), det, multiplies());
}

Matrix Matrix::inverse() const throw()
{
    assert(is_square());

    if( ! is_invertible())
        throw runtime_error("\aNOT INVERTIBLE\n");

    Matrix mx = *this;
    Matrix inverse = Matrix::Identity(rows_num);

    bool alternative_pivot_1_found;

    bool pivot_not_zero_found;

    bool number_not_zero_found;

    size_t row_with_alternative_pivot;

    size_t row_with_pivot_not_zero;

    size_t pivot_row = 0;
    size_t pivot_col = 0;

    //Gauss Elimination
    while (pivot_row < (rows_num - 1))
    {
        alternative_pivot_1_found = mx.pivotEqualTo_one_Found (pivot_row, pivot_col, row_with_alternative_pivot);

        pivot_not_zero_found = mx.pivotNot_zero_Found(pivot_row, pivot_col, row_with_pivot_not_zero);

        if (mx.at(pivot_row, pivot_col) != 1 && alternative_pivot_1_found )
        {
            inverse.swapRows(pivot_row, row_with_alternative_pivot);
            mx.swapRows(pivot_row, row_with_alternative_pivot);
        }
        else if (mx.at(pivot_row, pivot_col) == 0 && pivot_not_zero_found )
        {
            inverse.swapRows(pivot_row, row_with_pivot_not_zero);
            mx.swapRows(pivot_row, row_with_pivot_not_zero );
        }

        size_t col_dif_zero;

        number_not_zero_found = mx.firstNumberNot_zero(pivot_row, col_dif_zero);

        if(number_not_zero_found)
        {
            if ( mx.at(pivot_row, col_dif_zero) != 1)
            {
                inverse.changePivotTo_one(pivot_row, mx.at(pivot_row, col_dif_zero));
                mx.changePivotTo_one(pivot_row, mx.at(pivot_row, col_dif_zero));
            }
        }

        if(number_not_zero_found)
        {
            for (size_t i = pivot_row + 1; i < cols_num; ++i)
            {
                inverse.zeroOutTheColumn(i, pivot_row, mx.at(i, col_dif_zero));
                mx.zeroOutTheColumn(i, pivot_row, mx.at(i, col_dif_zero));
            }
        }

        ++pivot_row;
        ++pivot_col;
    }

    //Jordan Elimination
    while(pivot_row > 0)
    {
        size_t col_dif_zero;

        number_not_zero_found = mx.firstNumberNot_zero(pivot_row, col_dif_zero);

        if(number_not_zero_found)
        {
            if ( mx.at(pivot_row, col_dif_zero) != 1)
            {
                inverse.changePivotTo_one(pivot_row, mx.at(pivot_row, col_dif_zero));
                mx.changePivotTo_one(pivot_row, mx.at(pivot_row, col_dif_zero));
            }
        }

        if(number_not_zero_found)
        {
            for (int i = pivot_row - 1; i >= 0; --i)
            {
                inverse.zeroOutTheColumn(i, pivot_row, mx.at(i, col_dif_zero));
                mx.zeroOutTheColumn(i, pivot_row, mx.at(i, col_dif_zero));
            }
        }

        --pivot_row;
    }

    return inverse;
}

Matrix Matrix::adjoint() const
{
    assert(is_square());
    assert(cols_num > 1);

    if(is_zero())
        return Matrix(rows_num, cols_num);

    if(is_constant())
        return Matrix(rows_num, cols_num);

    if(is_identity())
        return *this;

    Matrix cofact(rows_num, cols_num);

    size_t r = 0, c = 0;

    Matrix temp(rows_num - 1, cols_num - 1);

    for(size_t i = 0; i < rows_num; ++i)
    {
        for(size_t j = 0; j < cols_num; ++j)
        {
            for(size_t k = 0; k < rows_num; ++k)
            {
                for(size_t h = 0; h < cols_num; ++h)
                {
                    if (k != i && h != j)
                    {
                        temp(r, c++) = at(k, h);

                        if(c == cols_num - 1)
                        {
                            c = 0;
                            ++r;
                        }
                    }
                }
            }

            c = 0;
            r = 0;

            int sign;

            sign = ( ( i + j ) % 2 == 0 ) ? 1 : -1;

            cofact.at(i, j) = sign * temp.determinant();
        }
    }

    return cofact.transpose();
}

Matrix Matrix::gaussJordanElimination() const
{
    Matrix mx = *this;

    bool alternative_pivot_1_found;

    bool pivot_not_zero_found;

    bool number_not_zero_found;

    size_t row_with_alternative_pivot;

    size_t row_with_pivot_not_zero;

    size_t pivot_row = 0;
    size_t pivot_col = 0;

    ///Gauss Elimination
    while (pivot_row < (rows_num - 1) && pivot_row < (cols_num))
    {
        alternative_pivot_1_found = mx.pivotEqualTo_one_Found ( pivot_row, pivot_col,
                                    row_with_alternative_pivot);

        pivot_not_zero_found = mx.pivotNot_zero_Found(
                                   pivot_row, pivot_col, row_with_pivot_not_zero);

        if (mx.at( pivot_row, pivot_col) != 1 && alternative_pivot_1_found )
        {
            mx.swapRows(pivot_row, row_with_alternative_pivot);
        }
        else if (mx.at( pivot_row, pivot_col) == 0 && pivot_not_zero_found )
        {
            mx.swapRows( pivot_row, row_with_pivot_not_zero );
        }

        size_t col_dif_zero;

        number_not_zero_found = mx.firstNumberNot_zero(pivot_row, col_dif_zero);

        if(number_not_zero_found)
        {
            if (( mx.at(pivot_row, col_dif_zero) ) != 1)
            {
                mx.changePivotTo_one(pivot_row,
                                     mx.at(pivot_row, col_dif_zero) );
            }
        }

        if(number_not_zero_found)
        {
            for(size_t i = pivot_row + 1; i < rows_num; ++i)
            {
                mx.zeroOutTheColumn( i, pivot_row, mx.at(i, col_dif_zero));
            }
        }

        ++pivot_row;
        ++pivot_col;
    }

    //Jordan Elimination
    while(pivot_row > 0)
    {
        size_t col_dif_zero;

        number_not_zero_found = mx.firstNumberNot_zero(pivot_row, col_dif_zero);

        if(number_not_zero_found)
            if ( mx.at(pivot_row, col_dif_zero) != 1)
            {
                mx.changePivotTo_one(pivot_row, mx.at(pivot_row, col_dif_zero));
            }

        if(number_not_zero_found)
            for (int i = pivot_row - 1; i >= 0; --i)
                mx.zeroOutTheColumn(i, pivot_row, mx.at(i, col_dif_zero));

        --pivot_row;
    }

    return mx;
}

Matrix Matrix::gaussElimination() const
{
    Matrix mx = *this;

    bool alternative_pivot_1_found;

    bool pivot_not_zero_found;

    bool number_not_zero_found;

    size_t row_with_alternative_pivot;

    size_t row_with_pivot_not_zero;

    size_t pivot_row = 0;
    size_t pivot_col = 0;

    ///Gauss Elimination
    while (pivot_row < (rows_num - 1) && pivot_row < (cols_num) )
    {
        alternative_pivot_1_found = mx.pivotEqualTo_one_Found ( pivot_row, pivot_col,
                                    row_with_alternative_pivot);

        pivot_not_zero_found = mx.pivotNot_zero_Found(
                                   pivot_row, pivot_col, row_with_pivot_not_zero);

        if (mx.at( pivot_row, pivot_col) != 1 && alternative_pivot_1_found )
        {
            mx.swapRows(pivot_row, row_with_alternative_pivot);
        }
        else if (mx.at( pivot_row, pivot_col) == 0 && pivot_not_zero_found )
        {
            mx.swapRows( pivot_row, row_with_pivot_not_zero );
        }

        size_t col_dif_zero;

        number_not_zero_found = mx.firstNumberNot_zero(pivot_row, col_dif_zero);

        if(number_not_zero_found)
        {
            if (( mx.at(pivot_row, col_dif_zero) ) != 1)
            {
                mx.changePivotTo_one(pivot_row,
                                     mx.at(pivot_row, col_dif_zero) );
            }
        }

        if(number_not_zero_found)
        {
            for(size_t i = pivot_row + 1; i < rows_num; ++i)
            {
                mx.zeroOutTheColumn( i, pivot_row, mx.at(i, col_dif_zero));
            }
        }

        ++pivot_row;
        ++pivot_col;
    }

    size_t col_dif_zero;

    number_not_zero_found = mx.firstNumberNot_zero(pivot_row, col_dif_zero);

    if(number_not_zero_found)
        if ( mx.at(pivot_row, col_dif_zero) != 1)
        {
            mx.changePivotTo_one(pivot_row, mx.at(pivot_row, col_dif_zero));
        }

    return mx;
}

vector<Fraction> Matrix::main_diagonal()
{
    assert(is_square());

    vector<Fraction> diag;

    for(size_t i = 0; i < rows_num; ++i)
        diag.push_back(at(i,i));

    return diag;
}

vector<Fraction> Matrix::secondary_diagonal()
{
    assert(is_square());

    vector<Fraction> diag;

    for(size_t i = 0, j = rows_num - 1; i < rows_num; ++i, --j)
        diag.push_back(at(i,j));

    return diag;
}

void Matrix::swapRows( int row1, int row2)
{
    for (size_t i = 0; i < cols_num; i++ )
        std::swap( at(row1,i ), at(row2, i) );
}

bool Matrix::pivotEqualTo_one_Found( int pivot_row, int pivot_col, size_t& alternative_pivot_row )
{
    for (size_t i = pivot_row + 1; i < rows_num; ++i)
    {
        if(at(i, pivot_col) == 1)
        {
            alternative_pivot_row = i;

            return true;
        }
    }

    return false;
}

bool Matrix::pivotNot_zero_Found(int pivot_row, int pivot_col,size_t& col_dif_zero )
{
    for (size_t i = pivot_row + 1; i < rows_num; ++i)
        if(at(i, pivot_col) != 0)
        {
            col_dif_zero = i;

            return true;
        }

    return false;
}

bool Matrix::firstNumberNot_zero(int row_num, size_t& num_coluna_num_dif_zero)
{
    for (size_t i = 0; i < cols_num; ++i)
        if (at(row_num, i) != 0)
        {
            num_coluna_num_dif_zero = i;

            return true;
        }

    return false;
}

void Matrix::changePivotTo_one( int row_num, Fraction constant)
{
    for(size_t i = 0; i < cols_num; ++i)
        if (at(row_num, i).num != 0)
            at(row_num, i) = (at(row_num, i) / constant);
}

void Matrix::zeroOutTheColumn( int row_num, int num_pivot_row, Fraction constant)
{
    for(size_t i = 0; i < cols_num; ++i)
        at(row_num, i) = at(row_num, i) -  (constant * at(num_pivot_row, i));
}

}// L_Algebra namespace
