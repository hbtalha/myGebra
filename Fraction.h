#ifndef FRACTION_H_INCLUDED
#define FRACTION_H_INCLUDED

#include <iostream>
#include <ostream>
#include <cstring>
#include <assert.h>
#include <stdint.h>

class Fraction
{
    using Integer = std::int64_t;

public:
    Integer num;
    Integer den;

    Fraction(Integer numm = 0, Integer denm = 1) : num(numm), den(denm)
    {
        if(den == 0)
            throw std::domain_error("denominator 0(zero) exception");
    }

    Fraction (Fraction n, Fraction d) : num(n.num * d.den), den(n.den * d.num)
    {
        if(den == 0)
            throw std::domain_error("denominator 0(zero) exception");
    }

    friend std::ostream& operator<< (std::ostream& os, const Fraction& fr);

    std::string to_string() const;

    bool is_integer()
    {
        return den == 1;
    }

    explicit operator bool() const
    {
        return num != 0;
    }

    bool operator== (const Fraction& fr) const
    {
        return num == fr.num && den == fr.den;
    }

    bool operator!= (const Fraction& fr) const
    {
        return !(*this == fr);
    }

    bool operator<(const Fraction& fr)
    {
        return num * fr.den < den * fr.num;
    }

    bool operator>(const Fraction& fr)
    {
        return !(*this <= fr);
    }

    bool operator<=(const Fraction& fr)
    {
        return *this == fr || *this < fr;
    }

    bool operator>=(const Fraction& fr)
    {
        return !(*this < fr);
    }

    Fraction operator-() const
    {
        return { -num, den };
    }

    Fraction operator+() const
    {
        return *this;
    }

    double to_double()
    {
        return double(num) / den;
    }

    float to_float()
    {
        return float(num) / den;
    }

    void simplify();

    friend Fraction operator+(const Fraction& lhs, const Fraction& rhs);
    friend Fraction operator-(const Fraction& lhs, const Fraction& rhs);
    friend Fraction operator*(const Fraction& lhs, const Fraction& rhs);
    friend Fraction operator/(const Fraction& lhs, const Fraction& rhs);

    void operator+= (const Fraction& fr);
    void operator-= (const Fraction& fr);
    void operator*= ( const Fraction& fr);
    void operator/= (const Fraction& fr);
};

Fraction pow_fract(const Fraction& fr, int x);

#endif // FRACTION_H_INCLUDED
