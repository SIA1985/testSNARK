#include "funciton.h"

snrk::Y_t snrk::Polynom::commit(TG_t tG)
{
    /*пока = f(t)*G напрямую*/
    return this->operator()(tG.t) * tG.G;
}

snrk::Lagrange snrk::Lagrange::generate(const dots_t &dots)
{
    Lagrange l;

    l.m_dots = dots;

    return l;
}

snrk::Y_t snrk::Lagrange::operator()(const X_t &x)
{
    Y_t y = 0;

    for(std::size_t i = 0; i < m_dots.size(); i++) {
        y += (l(i, x) * m_dots[i].y);
    }

    return y;
}

snrk::X_t snrk::Lagrange::l(std::size_t i, const X_t &x)
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

snrk::CustomPolynom snrk::CustomPolynom::generate(const func_t &customFunction)
{
    CustomPolynom c;

    c.m_customFunction = customFunction;

    return c;
}

snrk::Y_t snrk::CustomPolynom::operator()(const X_t &x)
{
    return m_customFunction(x);
}
