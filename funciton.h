#ifndef FUNCITON_H
#define FUNCITON_H

#include "circut.h"


namespace snrk {

template <typename X, typename Y>
struct DotValue {
    X x;
    Y y;
};

template <typename X, typename Y>
class polynom
{   
public:
    polynom() = default;
    ~polynom() = default;

    virtual Y operator()(const X &x) = 0;
};

template <typename X, typename Y>
class lagrange : public polynom<X, Y>
{
    using dot_t = DotValue<X, Y>;
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
        auto l = [this, &x](std::size_t i) -> X
        {
            X li;
            for(std::size_t j = 0; j < m_dots.size(); j++) {
                if (i == j) {
                    continue;
                }

                li *= ((x - m_dots[j].x) / (m_dots[i].x - m_dots[j].x));
            }

            return li;
        };

        Y y;

        for(std::size_t i = 0; i < m_dots.size(); i++) {
            y += (l(i) * m_dots[i].y);
        }

        return y;
    }

private:
    lagrange() = default;

    dots_t m_dots;
};

}
#endif // FUNCITON_H
