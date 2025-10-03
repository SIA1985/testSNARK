#ifndef FUNCITON_H
#define FUNCITON_H

#include "circut.h"

#include <functional>
#include <unordered_set>


namespace snrk {

template <typename X, typename Y>
class Polynom
{
public:
    Polynom() = default;
    virtual ~Polynom() = default;

    virtual Y operator()(const X &x) = 0;

    /*todo: передача gp*/
    virtual Y commit(X t, int G)
    {
        /*пока = f(t)*G напрямую*/
        return this->operator()(t) * G;

    }
};

template <typename X, typename Y>
class Lagrange : public Polynom<X, Y>
{
    using dot_t = struct{X x; Y y;};
    using dots_t = std::vector<dot_t>;

public:
    static Lagrange generate(const dots_t &dots)
    {
        Lagrange l;

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

private:
    Lagrange() = default;

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

template <typename X, typename Y>
class CustomPolynom : public Polynom<X, Y>
{
    using func_t = std::function<Y(X)>;
public:
    CustomPolynom(const func_t &customFunction)
        : m_customFunction{customFunction}
    {
    }

    virtual Y operator()(const X &x) override
    {
        return m_customFunction(x);
    }

private:
    func_t m_customFunction;
};

template <typename X, typename Y>
class ZeroPolynom : public Polynom<X, Y>
{
    using dot_t = struct{X x; Y y;};
    using dots_t = std::unordered_set<dot_t>;

public:
    static ZeroPolynom generate(const dots_t &dots)
    {
        ZeroPolynom z;

        z.m_dots = dots;

        return z;
    }

    virtual Y operator()(const X &x) override
    {
        return (m_dots.count(x) > 0 ? 0 : -1);
    }

private:
    ZeroPolynom() = default;

    dots_t m_dots;
};

}
#endif // FUNCITON_H
