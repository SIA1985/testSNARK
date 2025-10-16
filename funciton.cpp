#include "funciton.h"

#include <cmath>

namespace snrk {

Y_t Polynom::commit(TG_t tG)
{
    /*пока = f(t)*G напрямую*/
    return this->operator()(tG.t) * tG.G;
}

CustomPolynom CustomPolynom::generate(const func_t &customFunction)
{
    CustomPolynom c;

    c.m_customFunction = customFunction;

    return c;
}

Y_t CustomPolynom::operator()(const X_t &x)
{
    return m_customFunction(x);
}

ClassicPolynom ClassicPolynom::generate(coefs_t coefs)
{
    ClassicPolynom c;

    c.m_coefs = coefs;

    return c;
}

/*todo: тест*/
Y_t ClassicPolynom::operator()(const X_t &x)
{
    Y_t y = 0;

    for(std::size_t i = 0; i < m_coefs.size(); i++) {
        y += m_coefs[i] * std::pow(x, i);
    }

    return y;
}

/*todo: тест*/
CustomPolynom ClassicPolynom::operator/(ClassicPolynom &other)
{
    std::size_t k,j;

    std::size_t n = m_coefs.size() - 1;
    std::size_t nOther = other.m_coefs.size() - 1;

    ClassicPolynom res(n - nOther), rem(n - nOther);

    for(j = 0; j <= n; j++) {
            rem[j] = (*this)[j];
            res[j]=0.0;
    }
    for(k = n - nOther; k >= 0; k--) {
        res[k] = rem[nOther+k] / other[nOther];

        for(j = nOther + k - 1; j >= k; j--) {
            rem[j] -= res[k] * other[j-k];
        }
    }
    for(j = nOther; j <= n; j++) {
        rem[j] = 0.0;
    }

    return CustomPolynom::generate([res, rem](X_t x) mutable -> Y_t
    {
        return res(x) + rem(x);
    });
}

ValueType &ClassicPolynom::operator[](std::size_t i)
{
    return m_coefs[i];
}

ClassicPolynom::ClassicPolynom(std::size_t n)
{
    m_coefs.resize(n, ValueType(0));
}


Lagrange Lagrange::generate(const dots_t &dots)
{
    Lagrange l;

    l.m_dots = dots;

    return l;
}

Y_t Lagrange::operator()(const X_t &x)
{
    Y_t y = 0;

    for(std::size_t i = 0; i < m_dots.size(); i++) {
        y += (l(i, x) * m_dots[i].y);
    }

    return y;
}

X_t Lagrange::l(std::size_t i, const X_t &x)
{
    X_t li = 1;
    for(std::size_t j = 0; j < m_dots.size(); j++) {
        if (i == j) {
            continue;
        }

        li *= ((x - m_dots[j].x) / (m_dots[i].x - m_dots[j].x));
    }

    return li;
}

ZeroPolynom ZeroPolynom::generate(const xs_t &xs)
{
    ZeroPolynom z;

    for(const auto &x : xs) {
        z.m_dots.push_back({x, Y_t(0)});
    }

    return z;
}

}
