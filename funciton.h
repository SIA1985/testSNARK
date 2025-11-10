#ifndef FUNCITON_H
#define FUNCITON_H

#include "circut.h"

#include <functional>
#include <map>
#include <set>
#include <unordered_set>


namespace snrk {

typedef ValueType X_t;
typedef ValueType Y_t;

struct dot_t {X_t x; Y_t y;};
typedef std::vector<dot_t> dots_t;
struct TG_t{X_t t; int G;};
typedef std::set<X_t> xs_t;

class Polynom
{
public:
    using xProc_t = std::function<void(X_t&)>;

    Polynom() = default;
    virtual ~Polynom() = default;

    virtual Y_t operator()(X_t x) = 0;

    /*todo: передача gp*/
    virtual Y_t commit(TG_t tG);
};

class CustomPolynom : public Polynom
{
    using func_t = std::function<Y_t(X_t)>;
public:
    static CustomPolynom generate(const func_t &customFunction);

    virtual Y_t operator()(X_t x) override;

private:
    func_t m_customFunction;
};

class CanonicPolynom : public Polynom
{
    using coefs_t = std::vector<ValueType>;
    using roots_t = xs_t;
public:
    static CanonicPolynom generate(coefs_t coefs);

    static coefs_t coefsFromRoots(roots_t roots);

    virtual Y_t operator()(X_t x) override;

    CustomPolynom operator/(CanonicPolynom &other);

    CanonicPolynom operator+(const CanonicPolynom &other) const;

    CanonicPolynom operator-(CanonicPolynom &other);

    CanonicPolynom operator*(const CanonicPolynom &other) const;

    CanonicPolynom operator*(const ValueType value) const;

    void operator+=(const CanonicPolynom &other);

    void operator*=(const CanonicPolynom &other);

    CanonicPolynom operator()(const CanonicPolynom &other) const;

    ValueType &operator[](std::size_t i);

    std::size_t degree() const;

protected:
    CanonicPolynom() = default;
    CanonicPolynom(std::size_t n);

    /*x0, x1 .. xn*/
    coefs_t m_coefs;
};

/* [left, right] */
class Range
{
public:
    Range(X_t left, X_t right);

    int inRange(X_t x) const;

    X_t leftBound() const;
    X_t rightBound() const;

private:
    X_t m_left;
    X_t m_right;
};

bool operator<(const Range &a, const Range &b);

class PartedCanonicPloynom : public Polynom
{
    class RangeMap
    {
    public:
        CanonicPolynom operator[](X_t x);

        void insert(Range range, CanonicPolynom polynom);

    private:
        std::map<Range, CanonicPolynom> m_map;
    };

public:
    static PartedCanonicPloynom generate();

    virtual Y_t operator()(X_t x) override;

private:
    RangeMap m_map;
};

class InterpolationPolynom : public Polynom
{
    using coefs_t = std::vector<ValueType>;
public:
    static InterpolationPolynom generate(const dots_t &dots);

    virtual Y_t operator()(X_t x) override;

    CanonicPolynom toCanonicPolynom() const;
    PartedCanonicPloynom toPartedCanonicPolynom() const;

protected:
    dots_t m_dots;
};

class ZeroPolynom : public CanonicPolynom
{
public:
    static ZeroPolynom generate(const xs_t &xs);
};

}
#endif // FUNCITON_H
