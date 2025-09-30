#ifndef FUNCITON_H
#define FUNCITON_H

#include "circut.h"

#include <functional>


namespace snrk {

template <typename X, typename Y>
class polynom
{
public:
    polynom() = default;
    ~polynom() = default;

    virtual Y operator()(const X &x) = 0;

    /*todo: передача gp*/
    virtual Y commit(X t, int G) = 0;
};

template <typename X, typename Y>
class lagrange : public polynom<X, Y>
{
    using dot_t = struct{X x; Y y;};
    using dots_t = std::vector<dot_t>;

public:
    static lagrange generate(const dots_t &dots)
    {
        lagrange l;

        l.m_dots = dots;

        return l;
    }

    virtual Y operator()(const X &x) override
    {
        Y y = 0;

        for(std::size_t i = 0; i < m_dots.size(); i++) {
            y += (l(i, x) * m_dots[i].y);
        }

        return y;
    }

    virtual Y commit(X t, int G) override
    {
        /*пока = f(t)*G напрямую*/
        Y y = 0;
        for(std::size_t i = 0; i < m_dots.size(); i++) {
            y += (l(i, t) * G * m_dots[i].y);
        }

        return y;
    }

private:
    lagrange() = default;

    X l(std::size_t i, const X &x)
    {
        X li = 1;
        for(std::size_t j = 0; j < m_dots.size(); j++) {
            if (i == j) {
                continue;
            }

            li *= ((x - m_dots[j].x) / (m_dots[i].x - m_dots[j].x));
        }

        return li;
    };

    dots_t m_dots;
};

}
#endif // FUNCITON_H
