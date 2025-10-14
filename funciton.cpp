#include "funciton.h"

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

CustomPolynom polynomDevide(const Polynom &a, const Polynom &b)
{
    //1. При построении полинома найти корни (когда == 0)
    //2. Если корни совпадают, то удалить совпадения
    //3. Вернуть CustomPolynom
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

    z.m_xs = xs;

    return z;
}

Y_t ZeroPolynom::operator()(const X_t &x)
{
    return (m_xs.count(x) > 0 ? 0 : -1);
}

}
