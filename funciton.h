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


/*todo: порядок точек?*/
class Lagrange : public Polynom
{
public:
    static Lagrange generate(const dots_t &dots);

    virtual Y_t operator()(const X_t &x) override;

private:
    X_t l(std::size_t i, const X_t &x);

    dots_t m_dots;
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

//class ZeroPolynom : public Polynom
//{
//public:
//    static ZeroPolynom generate(const dots_t &dots)
//    {
//        ZeroPolynom z;

//        z.m_dots = dots;

//        return z;
//    }

//    virtual Y_t operator()(const X_t &x) override
//    {
//        return (m_dots.count(x) > 0 ? 0 : -1);
//    }

//private:
//    dots_t m_dots;
//};

}
#endif // FUNCITON_H
