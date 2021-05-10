#ifndef LA_VECTOR_H
#define LA_VECTOR_H

#include "Fraction.h"
#include "Matrix.h"
#include <initializer_list>
#include <deque>
#include <ostream>

namespace L_Algebra
{

class Vector
{
    std::vector<Fraction> data;

    Fraction& at(std::size_t i)
    {
        return data.at(i);
    }

    const Fraction& at(std::size_t i) const
    {
        return data.at(i);
    }

    void push_back(Fraction n)
    {
        data.push_back(n);
    }

    friend std::vector<Vector> null_space(Matrix mx);
    friend std::vector<Vector> null_space_(Matrix mx);

public:
    Vector() = default;

    template<std::input_iterator InputIt>
    Vector(InputIt first, InputIt last) : data(first, last) {}

    Vector(int d, long long int n = 0) : data(d, n) {}

    Vector(std::initializer_list<Fraction> values) : data(values) {}

    friend std::ostream& operator<< (std::ostream& os, const Vector& lav);

    explicit operator bool() const
    {
        return dimension() != 0;
    }

    bool operator==(const Vector& lav) const
    {
        return data == lav.data;
    }

    bool operator!=(const Vector& lav) const
    {
        return data != lav.data;
    }

    Fraction& operator[](size_t i)
    {
        return data.at(i);
    }

    const Fraction& operator[](size_t i) const
    {
        return data.at(i);
    }

    Vector operator+(const Vector& lav) const;
    Vector operator-(const Vector& lav) const;
    Vector operator->*(const Vector& lav) const; // vectorial product
    Fraction operator*(const Vector& lav) const; // dot product

    Vector& operator+=(const Vector& lav);
    Vector& operator-=(const Vector& lav);

    friend Vector operator*(const Vector& mx, Fraction n);
    friend Vector operator*(Fraction n, const Vector& mx);

    std::size_t dimension() const
    {
        return data.size();
    }

    Fraction norm_Power2() const;
    double norm() const;
};

Vector proj(Vector u, Vector a);
Vector proj_orthogonal(Vector u, Vector a);

bool is_orthogonal(std::initializer_list<Vector> vec_set);

bool is_linearly_dependent(std::initializer_list<Vector> vec_set);
bool is_linearly_dependent(std::initializer_list<Matrix> matrices_set);
bool is_linearly_independent(std::initializer_list<Vector> vec_set);
bool is_linearly_independent(std::initializer_list<Matrix> matrices_set);

bool is_linear_combination(std::initializer_list<Vector> vec_set, Vector vec);
bool is_linear_combination(std::initializer_list<Matrix> matrices_set, Matrix mx);

bool spans_space(std::initializer_list<Vector> vec_set);
bool spans_space(std::initializer_list<Matrix> matrix_set);
bool is_in_span(Vector vec, std::initializer_list<Vector> span);

bool is_basis(std::initializer_list<Vector> vec_set);
bool is_basis(std::initializer_list<Matrix> matrices_set);

Vector change_basis(Vector vec, std::initializer_list<Vector> basis_from, std::initializer_list<Vector> basis_to);
Vector change_basis(Vector vec_in_standard_basis, std::initializer_list<Vector> destination_basis);

std::vector<Vector> row_space_basis(Matrix mx);
std::vector<Vector> column_space_basis(Matrix mx);
std::vector<Vector> null_space(Matrix mx);
std::vector<Vector> null_space_(Matrix mx);

std::size_t row_space_dim(Matrix mx);
std::size_t column_space_dim(Matrix mx);
std::size_t nullity(Matrix mx);

Vector coordinate_vector_relative_to_basis(std::initializer_list<Vector> basis, Vector vec);
Vector vector_with_coordinate_relative_to_basis(std::initializer_list<Vector> basis, Vector coordinate_vec);

Matrix vectorsToMatrix(std::vector<Vector>vec_set);

Matrix turnMatricesIntoLinearCombination(std::vector<Matrix>matrix_set);

/*
Vector rowOfMatrixToVector(Matrix mx, int row);
Vector columnOfMatrixToVector(Matrix mx, int column);
*/

} // L_Algebra namespace

#endif // LA_VECTOR_H
