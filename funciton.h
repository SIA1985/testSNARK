#ifndef FUNCITON_H
#define FUNCITON_H

#include "circut.h"

#include <functional>
#include <set>
#include <unordered_set>


namespace snrk {

typedef ValueType X_t;
typedef ValueType Y_t;

struct dot_t {X_t x; Y_t y;};
typedef std::vector<dot_t> dots_t;
struct TG_t{X_t t; int G;};

class Polynom
{
public:
    Polynom() = default;
    virtual ~Polynom() = default;

    virtual Y_t operator()(const X_t &x) = 0;

    /*todo: передача gp*/
    virtual Y_t commit(TG_t tG);
};

class CustomPolynom : public Polynom
{
    using func_t = std::function<Y_t(X_t)>;
public:
    static CustomPolynom generate(const func_t &customFunction);

    virtual Y_t operator()(const X_t &x) override;

private:
    func_t m_customFunction;
};

CustomPolynom polynomDevide(const Polynom &a, const Polynom &b);

/*todo: порядок точек?*/
class Lagrange : public Polynom
{
public:
    static Lagrange generate(const dots_t &dots);

    virtual Y_t operator()(const X_t &x) override;

private:
    X_t l(std::size_t i, const X_t &x);

    dots_t m_dots;


    friend CustomPolynom polynomDevide(const Polynom &a, const Polynom &b);
};

class ZeroPolynom : public Polynom
{
    using xs_t = std::set<X_t>;
public:
    static ZeroPolynom generate(const xs_t &xs)
    {
        ZeroPolynom z;

        z.m_xs = xs;

        return z;
    }

    virtual Y_t operator()(const X_t &x) override
    {
        return (m_xs.count(x) > 0 ? 0 : -1);
    }

private:
    xs_t m_xs;
};

}
#endif // FUNCITON_H
