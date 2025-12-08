#ifndef POLYNOM_H
#define POLYNOM_H

#include "circut.h"

#include <functional>
#include <map>
#include <unordered_set>
#include <cassert>

#include "types.h"

namespace snrk {

bool operator<(const dot_t &a, const dot_t &b);


class Polynom
{
public:
    using xProc_t = std::function<void(X_t&)>;

    Polynom() = default;
    virtual ~Polynom() = default;

    virtual Y_t operator()(X_t x) = 0;
    dots_t dots(xs_t xs);

    /*todo: передача gp*/
    virtual Y_t commit(TG_t tG);
};

class CustomPolynom : public Polynom
{
public:
    using func_t = std::function<Y_t(X_t)>;


    CustomPolynom(const func_t &customFunction);

    virtual Y_t operator()(X_t x) override;

private:
    func_t m_customFunction;
};

class CanonicPolynom : public Polynom
{
public:
    using coefs_t = std::vector<ValueType>;
    using roots_t = xs_t;
    using devideResult_t = std::pair<CanonicPolynom, CanonicPolynom>;


    CanonicPolynom() = default;

    CanonicPolynom(coefs_t coefs);

    CanonicPolynom(std::size_t n);

    static coefs_t coefsFromRoots(roots_t roots);

    virtual Y_t operator()(X_t x) override;

    CustomPolynom operator/(CanonicPolynom &other);
    devideResult_t devide(CanonicPolynom &other);

    CanonicPolynom operator+(const CanonicPolynom &other) const;

    CanonicPolynom operator-(const CanonicPolynom &other) const;

    CanonicPolynom operator*(const CanonicPolynom &other) const;

    CanonicPolynom operator*(const ValueType value) const;

    void operator+=(const CanonicPolynom &other);

    void operator*=(const CanonicPolynom &other);

    CanonicPolynom operator()(const CanonicPolynom &other) const;

    ValueType &operator[](std::size_t i);
    const ValueType operator[](std::size_t i) const;

    std::size_t degree() const;

    bool isZero() const;

protected:
    static CanonicPolynom buildPolynomialRecursive(const roots_t& roots, roots_t::const_iterator start, roots_t::const_iterator end);

    /*x0, x1 .. xn*/
    coefs_t m_coefs;
};

/* [left, right] */
class Range
{
public:
    enum pos_t : int
    {
        left        = 1 << 0,
        right       = 1 << 1,
        crossed     = 1 << 2,
        equal       = 1 << 3,
        inside      = 1 << 4,
        outside     = 1 << 5,
    };

    Range(X_t left, X_t right);

    bool inRangeStrict(X_t x) const;

    bool inRange(X_t x) const;

    pos_t isCrossStrict(const Range &other) const;

    Range crossByStrict(const Range &other) const;

    X_t leftBound() const;
    X_t rightBound() const;

    static Range fromUnsorted(X_t a, X_t b);

private:
    X_t m_left;
    X_t m_right;
};

bool operator<(const Range &a, const Range &b);
bool operator==(const Range &a, const Range &b);
bool operator<=(const Range &a, const Range &b);
std::ostream &operator<<(std::ostream &out, const Range &r);

/*Непрерывный диапазон с одинаковым шагом*/
template<typename T>
class RangeMap
{
public:
    using const_iterator = typename std::map<Range, T>::const_iterator;

    RangeMap() = default;

    RangeMap(const_iterator begin, const_iterator end)
    {
        m_map = std::map<Range, T>(begin, end);
    }

    const_iterator find(X_t x) const
    {
        return std::lower_bound(m_map.begin(), m_map.end(), x,
        [](const std::pair<Range, T> &keyValue, X_t x) -> bool
        {
            return cmp(keyValue.first.rightBound(), x) == -1;
        });
    }

    T operator[](X_t x)
    {
        auto found = find(x);

        if (found != m_map.end()) {
            return found->second;
        }

        if (cmp(m_map.begin()->first.leftBound(), x) == 1) {
            return m_map.begin()->second;
        }

        return (--m_map.end())->second;
    }

    T &operator[](Range r)
    {
        if (r < m_map.begin()->first) {
            auto leftF = m_map.begin()->second;
            m_map[r] = leftF;
            return m_map[r];
        }

        if  ((--m_map.end())->first < r) {
            auto rigthF = (--m_map.end())->second;
            m_map[r] = rigthF;
            return m_map[r];
        }

        return m_map[r];
    }

    /*Непрерывное заполнение*/
    void insert(const Range &range, const T &polynom)
    {
        //todo:
//        std::cout <<  "cho " << range.leftBound() << " " <<  m_map.begin()->first.leftBound() << std::endl;
//        bool i1 = m_map.size() == 0;
//        bool i2 = range.leftBound() == (--m_map.end())->first.rightBound();
//        bool i3 = range.rightBound() == m_map.begin()->first.leftBound();
//        bool i4 = range == m_map.begin()->first;
//        bool i5 = range == (--m_map.end())->first;
//        assert(m_map.size() == 0 ||
//               range.leftBound() == (--m_map.end())->first.rightBound() ||
//               range.rightBound() == m_map.begin()->first.leftBound() ||
//               range == m_map.begin()->first ||
//               range == (--m_map.end())->first);

        m_map.insert({range, polynom});
    }

    const_iterator cbegin() const
    {
        return m_map.cbegin();
    }

    const_iterator cend() const
    {
        return m_map.cend();
    }

    std::size_t size() const
    {
        return m_map.size();
    }

private:
    std::map<Range, T> m_map;
};

class PartedCanonicPolynom : public Polynom
{
public:
    using map = RangeMap<CanonicPolynom>;


    PartedCanonicPolynom() = default;

    PartedCanonicPolynom(const RangeMap<CanonicPolynom> &map);

    PartedCanonicPolynom(const std::set<dot_t> &sortedDots, bool fromInterpolation = true);

    virtual Y_t operator()(X_t x) override;

    PartedCanonicPolynom operator+(const PartedCanonicPolynom &other) const;

    PartedCanonicPolynom operator-(const PartedCanonicPolynom &other) const;

    PartedCanonicPolynom operator*(const PartedCanonicPolynom &other) const;

    CustomPolynom operator/(PartedCanonicPolynom &other);

    PartedCanonicPolynom mustDevide(PartedCanonicPolynom &other) const;

    void operator+=(const PartedCanonicPolynom &other);

    Range distance() const;

    PartedCanonicPolynom cut(Range distance) const;

    const static int Partition;

protected:
    using operatorPred_t = std::function<void(RangeMap<CanonicPolynom>::const_iterator it,
                                              RangeMap<CanonicPolynom>::const_iterator itOther,
                                              Range currentRange)>;
    void operatorPrivate(const PartedCanonicPolynom &other, operatorPred_t pred) const;

     map m_map;
};

class InterpolationPolynom : public Polynom
{
public:
    using coefs_t = std::vector<ValueType>;


    InterpolationPolynom() = default;

    InterpolationPolynom(const dots_t &dots);

    virtual Y_t operator()(X_t x) override;

    CanonicPolynom toCanonicPolynom() const;
    PartedCanonicPolynom toPartedCanonicPolynom() const;

protected:
    dots_t m_dots;
};

class ZeroPolynom : public CanonicPolynom
{
public:
    using roots_t = xs_t;


    ZeroPolynom(const xs_t &xs);

    PartedCanonicPolynom toPartedCanonicPolynom() const;

private:
    roots_t m_roots;
};

}
#endif // POLYNOM_H
