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

class CanonicPolynom : public Polynom
{
    using coefs_t = std::vector<ValueType>;
public:
    static CanonicPolynom generate(coefs_t coefs);

    virtual Y_t operator()(const X_t &x) override;

    CustomPolynom operator/(CanonicPolynom &other);

    CanonicPolynom operator-(CanonicPolynom &other);

    ValueType &operator[](std::size_t i);

private:
    CanonicPolynom() = default;
    CanonicPolynom(std::size_t n);

    /*x0, x1 .. xn*/
    coefs_t m_coefs;
};


class InterpolationPolynom : public Polynom
{
    using coefs_t = std::vector<ValueType>;
public:
    static InterpolationPolynom generate(const dots_t &dots);

    virtual Y_t operator()(const X_t &x) override;

    CanonicPolynom toClassicPolynom() const;

protected:
    dots_t m_dots;
    coefs_t m_newtonCoefs;

    friend CustomPolynom polynomDevide(const Polynom &a, const Polynom &b);
};

class ZeroPolynom : public InterpolationPolynom
{
    using xs_t = std::set<X_t>;
public:
    static ZeroPolynom generate(const xs_t &xs);
};

}
#endif // FUNCITON_H
