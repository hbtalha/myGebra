#include "LA_Vector.h"

#include <iostream>
#include <math.h>
#include <assert.h>
#include <set>
#include <deque>
#include <algorithm>

using namespace std;

namespace L_Algebra
{

Matrix transitionMatrix(Matrix from, Matrix to)
{
    assert(from.size() == to.size());

    size_t rows_num = to.rows();
    size_t cols_num = to.cols();

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
        alternative_pivot_1_found = to.pivotEqualTo_one_Found (pivot_row, pivot_col, row_with_alternative_pivot);

        pivot_not_zero_found = to.pivotNot_zero_Found(pivot_row, pivot_col, row_with_pivot_not_zero);

        if (to.at(pivot_row, pivot_col) != 1 && alternative_pivot_1_found )
        {
            from.swapRows(pivot_row, row_with_alternative_pivot);
            to.swapRows(pivot_row, row_with_alternative_pivot);
        }
        else if (to.at(pivot_row, pivot_col) == 0 && pivot_not_zero_found )
        {
            from.swapRows(pivot_row, row_with_pivot_not_zero);
            to.swapRows(pivot_row, row_with_pivot_not_zero );
        }

        size_t col_dif_zero;

        number_not_zero_found = to.firstNumberNot_zero(pivot_row, col_dif_zero);

        if(number_not_zero_found)
        {
            if ( to.at(pivot_row, col_dif_zero) != 1)
            {
                from.changePivotTo_one(pivot_row, to.at(pivot_row, col_dif_zero));
                to.changePivotTo_one(pivot_row, to.at(pivot_row, col_dif_zero));
            }
        }

        if(number_not_zero_found)
        {
            for (int i = pivot_row + 1; i < cols_num; ++i)
            {
                from.zeroOutTheColumn(i, pivot_row, to.at(i, col_dif_zero));
                to.zeroOutTheColumn(i, pivot_row, to.at(i, col_dif_zero));
            }
        }

        ++pivot_row;
        ++pivot_col;
    }

    //Jordan Elimination
    while(pivot_row > 0)
    {
        size_t col_dif_zero;

        number_not_zero_found = to.firstNumberNot_zero(pivot_row, col_dif_zero);

        if(number_not_zero_found)
        {
            if ( to.at(pivot_row, col_dif_zero) != 1)
            {
                from.changePivotTo_one(pivot_row, to.at(pivot_row, col_dif_zero));
                to.changePivotTo_one(pivot_row, to.at(pivot_row, col_dif_zero));
            }
        }

        if(number_not_zero_found)
        {
            for (int i = pivot_row - 1; i >= 0; --i)
            {
                from.zeroOutTheColumn(i, pivot_row, to.at(i, col_dif_zero));
                to.zeroOutTheColumn(i, pivot_row, to.at(i, col_dif_zero));
            }
        }

        --pivot_row;
    }

    return from;
}

bool is_consistent(const Matrix& mx)
{
    int rows_num = mx.rows();
    int cols_num = mx.cols();

    Matrix mx1 = mx.gaussJordanElimination();

    bool square = mx.is_square();

    int num_non_zero_numbers = 0;
    for(int i = 0; i < rows_num; ++i)
    {
        if (square)
            for(int j = 0; j < cols_num; ++j)
            {
                if(mx1(i, j) != 0)
                    ++num_non_zero_numbers;
            }
        else
            for(int j = 0; j < cols_num - 1; ++j)
            {
                if(mx1(i, j) != 0)
                    ++num_non_zero_numbers;
            }

        if(num_non_zero_numbers > 1)
            return false;

        if( ! square && num_non_zero_numbers == 0 && mx1(i, cols_num - 1) != 0)
            return false;

        num_non_zero_numbers = 0;
    }

    return true;
}

Matrix vectorsToMatrix(std::vector<Vector>vec_set)
{
    assert(vec_set.size() > 0);

    int len = vec_set.size();
    for(int i = 0; i < len; ++i)
        assert(vec_set[i].dimension() == vec_set[0].dimension());

    int rows_num = vec_set[0].dimension();
    int cols_num = len;

    Matrix mx(rows_num, cols_num);

    for(int i = 0; i < rows_num; ++i)
        for(int j = 0; j < cols_num; ++j)
        {
            mx(i, j) = vec_set.at(j)[i];
        }

    return mx;
}

Matrix turnMatricesIntoLinearCombination(std::vector<Matrix>matrix_set)
{
    assert(matrix_set.size() > 0);

    int len = matrix_set.size();
    for(int i = 0; i < len; ++i)
        assert(matrix_set[i].size() == matrix_set[0].size());
    /*
        int rows_num = matrix_set[0].size();
        int cols_num = len;

        int r = matrix_set[0].rows();
        int c = matrix_set[0].cols();

        Matrix m(rows_num, cols_num);

        Vector lav(r * c);

        size_t vec_lav_size = cols_num;
        vector<Vector> vec_lav(vec_lav_size, r * c);

        // pass the values from the set of matrices to a set of la_vectors
        int ind = 0;
        for(size_t h = 0; h < vec_lav_size; ++h)
        {
            for(int i = 0; i < r; ++i)
                for(int j = 0; j < c; ++j)
                    vec_lav.at(h)[ind++] = matrix_set.at(h)(i, j);

            ind = 0;
        }

         transform the values from the set of the matrices into a new matrix;
        for(int i = 0; i < rows_num; ++i)
            for(int j = 0; j < cols_num; ++j)
                m(i, j) = vec_lav.at(j)[i];
    */

    int rows_num = matrix_set[0].size();
    int cols_num = len;

    int r = matrix_set[0].rows();
    int c = matrix_set[0].cols();

    Matrix mtrx(rows_num, cols_num);

    for(int i = 0; i < cols_num; ++i)
    {
        int id = 0;

        for(int x = 0; x < r; ++x)
        {
            for(int y = 0; y < c; ++y)
            {
                mtrx(id++, i) = matrix_set[ i ](x, y);
            }
        }
    }

    return mtrx;
}

Vector rowOfMatrixToVector(const Matrix& mx, int row)
{
    assert(row <= mx.rows());

    int cols_num = mx.cols();

    Vector v(cols_num);

    for(int i = 0; i < cols_num; ++i)
        v[ i ] = mx(row, i);

    return v;
}

Vector columnOfMatrixToVector(const Matrix& mx, int column)
{
    assert(column <= mx.cols());

    int rows_num = mx.rows();

    Vector v(rows_num);

    for(int i = 0; i < rows_num; ++i)
        v[ i ] = mx(i, column);

    return v;
}

ostream& operator<< (ostream& os, const Vector& lav)
{
    os << "(";

    for(auto el : lav.data)
        os << el << ", ";

    if(lav.data.empty())
        os << " )";
    else
        os << "\b\b \b" << ")";

    return os;
}

Vector Vector::operator+(const Vector& lav) const
{
    size_t len = data.size();

    assert(len == lav.data.size());

    Vector addition;

    addition.data.resize(len, 0);

    for(size_t i = 0; i < len; ++i)
        addition[i] = at(i) + lav[i];

    return addition;
}

Vector& Vector::operator+=(const Vector& lav)
{
    return *this = *this + lav;
}

Vector Vector::operator-(const Vector& lav) const
{
    size_t len = data.size();

    assert(len == lav.data.size());

    Vector subtraction;

    subtraction.data.resize(data.size(), 0);

    for(size_t i = 0; i < len; ++i)
        subtraction[i] = at(i) - lav[i];

    return subtraction;
}

Vector& Vector::operator-=(const Vector& lav)
{
    return *this = *this - lav;
}

Fraction Vector::operator*(const Vector& lav) const // dot product
{
    size_t len = data.size();

    assert(len == lav.data.size());

    Fraction dot_prod;

    for(size_t i = 0; i < len; ++i)
        dot_prod += at(i) * lav[i];

    return dot_prod;
}

// vectorial product
Vector Vector::operator->*(const Vector& lav) const
{
    size_t len = data.size();

    assert( (len == lav.data.size()) && len == 3);

    return {at(1) * lav.at(2) - at(2) * lav.at(1),
            - (at(2) * lav.at(0) - at(0) * lav.at(2)),
            at(0) * lav.at(1) - at(1) * lav.at(0) };
}

Vector operator*(const Vector& lav, Fraction n)
{
    Vector mult;

    mult.data.resize(lav.data.size(), 0);

    int i = 0;
    for( auto el : lav.data)
        mult.at(i++) = el * n;

    return mult;
}

Vector operator*(Fraction n, const Vector& lav)
{
    Vector mult;

    mult.data.resize(lav.data.size(), 0);

    int i = 0;
    for( auto el : lav.data)
        mult.at(i++) = el * n;

    return mult;
}

double Vector::norm() const
{
    Fraction n;

    size_t len = dimension();

    for(size_t i = 0; i < len; ++i)
        n += pow_fract(at(i), 2);

    return sqrt(n.to_double());
}

Fraction Vector::norm_Power2() const
{
    Fraction n;

    size_t len = dimension();

    for(size_t i = 0; i < len; ++i)
        n += pow_fract(at(i), 2);

    return n;
}

bool is_orthogonal(std::initializer_list<Vector> vec_set)
{
    assert(vec_set.size() > 1);

    std::vector<Vector> vec(vec_set);

    size_t len = vec.size();

    for(size_t i = 0; i < len; ++i )
        assert(vec.at(i).dimension() == vec.at(0).dimension());

    for( size_t i = 0; i < len - 1; ++i)
        for( size_t j = i + 1; j < len; ++j)
            if (vec.at(i) * vec.at(j) == 0)
                return true;

    return false;
}

Vector proj(Vector u, Vector a)
{
    return Fraction(u*a, a.norm_Power2()) * a;
}

Vector proj_orthogonal(Vector u, Vector a)
{
    return u - proj(u, a);
}

bool is_linearly_dependent(std::initializer_list<Vector> vec_set)
{
    Matrix mx = vectorsToMatrix(vec_set).gaussJordanElimination();

    int rows_num = mx.rows();
    int cols_num = mx.cols();

    int num_non_zero_numbers = 0;
    for(int i = 0; i < rows_num; ++i)
    {
        for(int j = 0; j < cols_num; ++j)
        {
            if(mx(i, j) != 0)
                ++num_non_zero_numbers;
        }

        if(num_non_zero_numbers > 1)
            return true;

        num_non_zero_numbers = 0;
    }

    return false;
}

bool is_linearly_dependent(initializer_list<Matrix> matrices_set)
{
    assert(matrices_set.size() > 0);

    vector<Matrix> vecs(matrices_set);

    int len = vecs.size();
    for(int i = 0; i < len; ++i)
        assert(vecs[i].size() == vecs[0].size() && vecs[i].size() > 0);

    int r = vecs[0].rows();
    int c = vecs[0].cols();

    Matrix mx(r, c);

    vecs.push_back(mx);

    Matrix m = turnMatricesIntoLinearCombination(vecs);

    if( is_consistent(m))
        return false;
    else
        return true;
}

bool is_linearly_independent(std::initializer_list<Vector>vec_set)
{
    return ! is_linearly_dependent(vec_set);
}

bool is_linearly_independent(initializer_list<Matrix> matrices_set)
{
    return ! is_linearly_dependent(matrices_set);
}

bool is_linear_combination(std::initializer_list<Vector> vec_set, Vector vec)
{
    vector<Vector> vecs(vec_set);

    vecs.push_back(vec);

    Matrix mx = vectorsToMatrix(vecs);

    if( ! is_consistent(mx))
        return false;

    mx = mx.gaussJordanElimination();

    int rows_num = mx.rows();
    int cols_num = mx.cols();

    Vector results = columnOfMatrixToVector(mx, cols_num - 1);

    Vector combination(rows_num);

    for(int i = 0; i < rows_num; ++i)
    {
        for(int j = 0; j < cols_num - 1; ++j)
            combination[i] += results[j] * vecs.at(j)[i];
    }

    if(vec == combination)
        return true;
    else
        return false;
}

bool is_linear_combination(std::initializer_list<Matrix> matrices_set, Matrix mx)
{
    assert(matrices_set.size() > 0);

    vector<Matrix> vecs(matrices_set);
    vecs.push_back(mx);

    Matrix m = turnMatricesIntoLinearCombination(vecs);

    int cols_num = m.cols();

    vector<Vector> vec_lav(cols_num);

    for(int i = 0; i < cols_num; ++i)
        vec_lav[i] = columnOfMatrixToVector(m, i);

    if( ! is_consistent(m))
        return false;

    m = m.gaussJordanElimination();

    Vector results = columnOfMatrixToVector(m, cols_num - 1);

    Vector combination(m.rows());

    for(int i = 0; i < cols_num - 1; ++i)
        combination += results[i] * vec_lav.at(i);

    Vector lav = vec_lav[vec_lav.size() - 1];

    if(lav == combination)
        return true;
    else
        return false;
}

bool is_basis(std::initializer_list<Vector> vec_set)
{
    assert(vec_set.size() > 0);

    vector<Vector> vec(vec_set);

    int len = vec.size();
    for(int i = 0; i < len; ++i)
        assert(vec[i].dimension() == vec[0].dimension());

    if(vec.size() != vec[0].dimension())
        return false;

    return ! is_linearly_dependent(vec_set);
}

bool is_basis(std::initializer_list<Matrix> matrices_set)
{
    return ! is_linearly_dependent(matrices_set);
}

Vector change_basis(Vector vec, std::initializer_list<Vector> basis_from,
                    std::initializer_list<Vector> basis_to)
{
    assert(basis_to.size() == basis_from.size());
    assert(vec.dimension() == basis_from.size());

    Matrix from = vectorsToMatrix(basis_from);
    Matrix to = vectorsToMatrix(basis_to);

    Matrix transition_matrix = transitionMatrix(from, to);

    int vec_dimension = vec.dimension();

    Matrix vec_matrix(vec_dimension, 1);

    for(int i = 0; i < vec_dimension; ++i)
        vec_matrix(i,0) = vec[i];

    Matrix new_basis_vec_matrix = transition_matrix * vec_matrix;

    Vector vec_in_new_basis(vec_dimension);

    for(int i = 0; i < vec_dimension; ++i)
        vec_in_new_basis[i] = new_basis_vec_matrix(i,0);

    return vec_in_new_basis;
}

Vector change_basis(Vector vec_in_standard_basis, std::initializer_list<Vector> destination_basis)
{
    return coordinate_vector_relative_to_basis(destination_basis, vec_in_standard_basis);
}

bool spans_space(std::initializer_list<Vector> vec_set)
{
    return ! is_linearly_dependent(vec_set);
}

bool spans_space(std::initializer_list<Matrix> matrix_set)
{
    return ! is_linearly_dependent(matrix_set);
}

bool is_in_span(Vector vec, std::initializer_list<Vector> span)
{
    return is_linear_combination(span, vec);
}

Vector coordinate_vector_relative_to_basis(std::initializer_list<Vector> basis,
        Vector vec)
{
    assert(basis.size() == vec.dimension());

    vector<Vector> vecs(basis);

    vecs.push_back(vec);

    Matrix mx = vectorsToMatrix(vecs);

    mx = mx.gaussJordanElimination();

    int rows_num = mx.rows();
    int cols_num = mx.cols();

    if(! is_consistent(mx))
        throw runtime_error("the basis is linearly dependent");

    Vector coordinate_vector(rows_num);

    for(int i = 0; i < rows_num; ++i)
        coordinate_vector[i] = mx(i, cols_num - 1);

    return coordinate_vector;
}

Vector vector_with_coordinate_relative_to_basis(initializer_list<Vector> basis,
        Vector coordinate_vec)
{
    assert(basis.size() > 0);

    assert(coordinate_vec.dimension() == basis.size());

    vector<Vector> vecs(basis);

    int len = vecs.size();
    for(int i = 0; i < len; ++i)
        assert(vecs[i].dimension() == vecs[0].dimension());

    assert(coordinate_vec.dimension() == vecs[0].dimension());

    size_t basis_size = basis.size();
    size_t vec_size = vecs[0].dimension();

    Vector vec(vec_size);

    for(size_t i = 0; i < basis_size; ++i)
        for(size_t j = 0; j < vec_size; ++j)
            vec[i] += coordinate_vec[j] * vecs.at(j)[i];

    return vec;
}

std::vector<Vector> row_space_basis(Matrix mx)
{
    mx = mx.gaussJordanElimination();

    int rows_num = mx.rows();
    int cols_num = mx.cols();

    vector<Vector> space_basis;
    Vector lav(cols_num);

    for(int i = 0; i < rows_num; ++i)
        for(int j = 0; j < cols_num; ++j)
            if(mx(i, j) != 0)
            {
                for(int j = 0; j < cols_num; ++j)
                    lav[j] = mx(i, j);

                space_basis.push_back(lav);

                break;
            }

    return space_basis;
}

vector<Vector> column_space_basis(Matrix mx)
{
    Matrix m = mx.gaussJordanElimination();

    int rows_num = mx.rows();
    int cols_num = mx.cols();

    vector<Vector> space_basis;

    for(int i = 0; i < rows_num; ++i)
        for(int j = 0; j < cols_num; ++j)
        {
            Vector temp(rows_num);

            if(m(i, j) != 0)
            {
                for(int k = 0; k < rows_num; ++k)
                    temp[ k ] = mx(k, j);

                space_basis.push_back(temp);

                break;
            }
        }

    return space_basis;
}

vector<Vector> null_space_(Matrix mx)
{
    Matrix m = mx.gaussJordanElimination();

    int rows_num = mx.rows();
    int cols_num = mx.cols();

    vector<int> pivot_cols;

    cout << m << endl;

    for(int i = 0; i < rows_num; ++i)
        for(int j = 0; j < cols_num; ++j)
            if(m(i, j) != 0)
            {
                // keeps all cols numbers so it is guaranteed that the column that contains a pivot won't
                // be used for the null space
                pivot_cols.push_back(j);

                break;
            }

    vector<Vector> free_variables;//(cols_num);

    for(int i = 0; i < cols_num; ++i)
    {
        if(find(pivot_cols.begin(), pivot_cols.end(), i) == pivot_cols.end())
        {
            Vector lvec;
            for(int j = 0; j < rows_num; ++j)
                if(find(pivot_cols.begin(), pivot_cols.end(), i) == pivot_cols.end())
                    lvec.push_back(- m(j, i));


            if(find(pivot_cols.begin(), pivot_cols.end(), i) == pivot_cols.end())
                free_variables.push_back(lvec);
        }
    }

    int ind = 0;
    for(int i = rows_num; i < cols_num; ++i)
    {
        for(size_t j = 0; j < free_variables.size(); ++j)
            if(j == ind)
                free_variables[j].push_back(1);
            else
                free_variables[j].push_back(0);

        ++ind;
    }

    return free_variables;
}

vector<Vector> null_space(Matrix mx)
{
    Matrix m = mx.gaussJordanElimination();

    int rows_num = mx.rows();
    int cols_num = mx.cols();

    vector<int> pivot_cols;

    vector<Vector> free_variables(cols_num);

    for(int i = 0; i < rows_num; ++i)
        for(int j = 0; j < cols_num; ++j)
            if(m(i, j) != 0)
            {
                // keeps all cols numbers so it is guaranteed that the column that contains a pivot won't
                // be used for the null space
                pivot_cols.push_back(j);

                break;
            }

    for(int i = 0; i < rows_num; ++i)
        for(int j = 0; j < cols_num; ++j)
        {
            if(m(i,j) != 0)
            {
                for(int k = 0; k < cols_num; ++k)
                {
                    // the j'th column is the one with pivot so it can not be used for the null space
                    // meaning that it has to be above or below

                    // if it is below it means that the k'th column might be one with free variable,
                    // it will be checked, if it is free it will be added zero because to get to the
                    // j'th column it had to get past only zeroes
                    if( k < j )
                    {
                        // starting from the second row, before immediately adding 0(zero), it will be checked
                        // whether the column is one that contains a pivot, in case it does the 0 won't be added
                        if(i > 0)
                        {
                            if(find(pivot_cols.cbegin(), pivot_cols.cend(), k) == pivot_cols.cend())
                                free_variables[j].push_back(0);
                        }
                        else
                            free_variables[j].push_back(0);
                    }
                    else if(k > j && find(pivot_cols.cbegin(), pivot_cols.cend(), k) == pivot_cols.cend())
                    {
                        free_variables[j].push_back( -m(i, k) );
                    }
                }
                break;
            }
        }

    int num_vectors = free_variables.size();
    int dimension;

    // get the dimension of the vector that will be of the null space
    for(int i = 0; i < num_vectors; ++i)
        if (free_variables[i].dimension() != 0)
        {
            dimension = free_variables[i].dimension();
            break;
        }


    // add the Identity Matrix to the rows in the new matrix which correspond to the 'free' columns
    // in the original matrix, making sure the number of rows equals the number of columns in the
    // original matrix (otherwise, we couldn't multiply the original matrix against our new matrix)
    int ind = 0;
    for(int i = 0; i < num_vectors; ++i)
    {
        if(free_variables[i].dimension() == 0)
        {
            for(int j = 0; j < dimension; ++j)
                if(j == ind)
                    free_variables[i].push_back(1);
                else
                    free_variables[i].push_back(0);

            ++ind;
        }
    }

    vector<Vector> space_basis(dimension, num_vectors);

    for(int i = 0; i < dimension; ++i)
        for(int j = 0; j < num_vectors; ++j)
            space_basis.at(i)[ j ] = free_variables.at(j)[i];

    return space_basis;
}

std::size_t column_space_dim(Matrix mx)
{
    mx = mx.gaussJordanElimination();

    int rows_num = mx.rows();
    int cols_num = mx.cols();
    int dimension = 0;

    for(int i = 0; i < rows_num; ++i)
        for(int j = 0; j < cols_num; ++j)
            if(mx(i, j) != 0)
            {
                ++dimension;

                break;
            }

    return dimension;
}

std::size_t row_space_dim(Matrix mx)
{
    mx = mx.gaussJordanElimination();

    int rows_num = mx.rows();
    int cols_num = mx.cols();
    int dimension = 0;

    for(int i = 0; i < rows_num; ++i)
        for(int j = 0; j < cols_num; ++j)
            if(mx(i, j) != 0)
            {
                ++dimension;
                break;
            }

    return dimension;
}

std::size_t nullity(Matrix mx)
{
    Matrix m = mx.gaussJordanElimination();

    int rows_num = mx.rows();
    int cols_num = mx.cols();

    vector<int> pivot_cols;

    vector<Vector> free_variables(cols_num);

    for(int i = 0; i < rows_num; ++i)
        for(int j = 0; j < cols_num; ++j)
            if(m(i, j) != 0)
            {
                pivot_cols.push_back(j);

                break;
            }

    int dimension = 0;

    for(int i = 0; i < rows_num; ++i)
        for(int j = 0; j < cols_num; ++j)
            if(m(i,j) != 0)
            {
                for(int k = 0; k < cols_num; ++k)
                {
                    if(k < j )
                    {
                        if(i > 0)
                        {
                            if(find(pivot_cols.cbegin(), pivot_cols.cend(), k) == pivot_cols.cend())
                                ++dimension;
                        }
                        else
                            ++dimension;
                    }
                    else if(k > j && find(pivot_cols.cbegin(), pivot_cols.cend(), k) == pivot_cols.cend())
                        ++dimension;
                }

                return dimension;
            }

    return 0;
}

}// L_Algebra namespace
