#include "Fraction.h"

#include <iostream>
#include <ostream>
#include <sstream>

using namespace std;

using Integer = int64_t;

std::ostream& operator<< (std::ostream& os, const Fraction& fr)
{
    if(fr.den == 1)
        os << fr.num;
    else
        os << fr.num << "/" << fr.den;

    return os;
}

string Fraction::to_string() const
{
    ostringstream os;
    os << *this;
    return os.str();
}

Integer gcf(Integer a, Integer b)
{
    if( b == 0)
        return abs(a);
    else
        return gcf(b, a%b);
}

void Fraction::simplify()
{
    if (den == 0 || num == 0)
    {
        num = 0;
        den = 1;
    }
    // Put neg. sign in numerator only.
    if (den < 0)
    {
        num *= -1;
        den *= -1;
    }

    // Factor out GCF from numerator and denominator.
    Integer n = gcf(num, den);
    num = num / n;
    den = den / n;
}

Fraction operator-(const Fraction& lhs, const Fraction& rhs)
{
    Fraction fr((lhs.num * rhs.den) - (rhs.num * lhs.den), lhs.den * rhs.den );

    fr.simplify();

    return fr;
}

Fraction operator+(const Fraction& lhs, const Fraction& rhs)
{
    Fraction fr((lhs.num * rhs.den) + (rhs.num * lhs.den), lhs.den * rhs.den );

    fr.simplify();

    return fr;
}

Fraction operator*(const Fraction& lhs, const Fraction& rhs)
{
    Fraction fr(lhs.num * rhs.num, lhs.den * rhs.den);

    fr.simplify();

    return fr;
}

Fraction operator/(const Fraction& lhs, const Fraction& rhs)
{
    Fraction fr (lhs.num * rhs.den, lhs.den * rhs.num);

    fr.simplify();

    return fr;
}

void Fraction::operator+=(const Fraction& fr)
{
    *this = *this + fr;
}

void Fraction::operator-=(const Fraction& fr)
{
    *this = *this - fr;
}

void Fraction::operator/=(const Fraction& fr)
{
    *this = *this / fr;
}

void Fraction::operator*=(const Fraction& fr)
{
    *this = *this * fr;
}

Fraction pow_fract(const Fraction& fr, int x)
{
    Fraction p(fr);

    for(int i = 0; i < x - 1; ++i)
        p *= fr;

    return p;
}
