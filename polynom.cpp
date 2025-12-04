#include "polynom.h"

#include <cmath>
#include <iostream>

namespace snrk {

bool operator<(const dot_t &a, const dot_t &b)
{
    return cmp(a.x, b.x) != 1;
}

dots_t Polynom::dots(xs_t xs)
{
    dots_t result;
    for(const auto &x : xs) {
        result.push_back({x, (*this)(x)});
    }

    return result;
}

Y_t Polynom::commit(TG_t tG)
{
    /*пока = f(t)*G напрямую*/
    return this->operator()(tG.t) * tG.G;
}

CustomPolynom::CustomPolynom(const func_t &customFunction)
    : m_customFunction{customFunction}
{
}

Y_t CustomPolynom::operator()(X_t x)
{
    return m_customFunction(x);
}

CanonicPolynom::CanonicPolynom(coefs_t coefs)
    : m_coefs{coefs}
{
}

CanonicPolynom::CanonicPolynom(std::size_t n)
{
    m_coefs.resize(n, ValueType(0.0));
}

CanonicPolynom::coefs_t CanonicPolynom::coefsFromRoots(roots_t roots)
{
    CanonicPolynom result(1);
    result[0] = 1.0;

    for (const auto &root : roots) {
        CanonicPolynom binomial(2);
        binomial[0] = -root;
        binomial[1] = 1.0;

        result = result * binomial;
    }

    return result.m_coefs;
}

Y_t CanonicPolynom::operator()(X_t x)
{
    Y_t y = 0;
    X_t xPow = 1;

    for(std::size_t i = 0; i < m_coefs.size(); i++) {
        y += xPow * m_coefs[i];
        xPow *= x;
    }

    return y;
}

CustomPolynom CanonicPolynom::operator/(CanonicPolynom &other)
{
    int k,j;

    std::size_t n = m_coefs.size() - 1;
    std::size_t nOther = other.m_coefs.size() - 1;

    /*todo: deg(res) = n - nOther, deg(rem) = nOther - 1*/
    CanonicPolynom res(n + 1), rem(n + 1);

    if (nOther < 0 || (nOther == 0 && std::abs(other.m_coefs[0].get_d()) < 1e-9)) {
        throw std::invalid_argument("Division by zero polynomial or zero constant.");
    }

    if (n < nOther) {
        res = CanonicPolynom(1);
        res[0] = 0.0;
        rem = *this;
        return CustomPolynom([rem, other](X_t x) mutable -> Y_t
        {
            return (rem(x) == 0 ? Y_t(0) : rem(x) / other(x));
        });
    }

    coefs_t remainderCoefs = m_coefs;

    int degreeRes = n - nOther;
    res = CanonicPolynom(degreeRes + 1);
    for(int i = 0; i <= degreeRes; ++i) {
        res[i] = 0.0;
    }

    for (int i = degreeRes; i >= 0; --i) {
        ValueType factor = remainderCoefs[i + nOther] / other.m_coefs[nOther];
        res[i] = factor;

        for (int j = 0; j <= nOther; ++j) {
            remainderCoefs[i + j] -= factor * other.m_coefs[j];
        }
    }

    rem = CanonicPolynom(nOther);
    for(int j = 0; j < nOther; ++j) {
         rem[j] = remainderCoefs[j];
    }

    return CustomPolynom([res, rem, other](X_t x) mutable -> Y_t
    {
        return res(x) + (rem(x) == 0 ? Y_t(0) : rem(x) / other(x));
    });
}

CanonicPolynom CanonicPolynom::operator+(const CanonicPolynom &other) const
{
    coefs_t diff;
    size_t max_size = std::max(m_coefs.size(), other.m_coefs.size());
    diff.resize(max_size, 0.0);

    for (size_t i = 0; i < max_size; ++i) {
        ValueType c1 = (i < m_coefs.size()) ? m_coefs[i] : 0.0;
        ValueType c2 = (i < other.m_coefs.size()) ? other.m_coefs[i] : 0.0;
        diff[i] = c1 + c2;
    }

    return CanonicPolynom(diff);
}

CanonicPolynom CanonicPolynom::operator-(const CanonicPolynom &other) const
{
    coefs_t diff;
    size_t max_size = std::max(m_coefs.size(), other.m_coefs.size());
    diff.resize(max_size, 0.0);

    for (size_t i = 0; i < max_size; ++i) {
        ValueType c1 = (i < m_coefs.size()) ? m_coefs[i] : 0.0;
        ValueType c2 = (i < other.m_coefs.size()) ? other.m_coefs[i] : 0.0;
        diff[i] = c1 - c2;
    }

    return CanonicPolynom(diff);
}

CanonicPolynom CanonicPolynom::operator*(const CanonicPolynom &other) const
{
    std::size_t newDegree = m_coefs.size() + other.m_coefs.size() - 1;
    CanonicPolynom result(newDegree);

    for(std::size_t i = 0; i < m_coefs.size(); i++) {
        for(std::size_t j = 0; j < other.m_coefs.size(); j++) {
            result[i + j] += m_coefs[i] * other.m_coefs[j];
        }
    }
    return result;
}

CanonicPolynom CanonicPolynom::operator*(const ValueType value) const
{
    coefs_t result(m_coefs.size(), 0);

    for(std::size_t i = 0; i < m_coefs.size(); i++) {
        result[i] = m_coefs[i] * value;
    }

    return CanonicPolynom(result);
}

void CanonicPolynom::operator+=(const CanonicPolynom &other)
{
    m_coefs = (*this + other).m_coefs;
}

void CanonicPolynom::operator*=(const CanonicPolynom &other)
{
    m_coefs = (*this * other).m_coefs;
}

CanonicPolynom CanonicPolynom::operator()(const CanonicPolynom &other) const
{
    if (m_coefs.size() <= 0) {
        return CanonicPolynom({m_coefs.back()});
    }

    CanonicPolynom result = CanonicPolynom({m_coefs.back()});

    for (int i = m_coefs.size() - 2; i >= 0; --i) {
        result *= other;
        result += CanonicPolynom({m_coefs[i]});
    }

    return result;
}

ValueType &CanonicPolynom::operator[](std::size_t i)
{
    return m_coefs[i];
}

const ValueType CanonicPolynom::operator[](std::size_t i) const
{
    return m_coefs.size() < i ? 0.0 : m_coefs.at(i);
}

std::size_t CanonicPolynom::degree() const
{
    return m_coefs.size() - 1;
}

InterpolationPolynom::InterpolationPolynom(const dots_t &dots)
    : m_dots{dots}
{
}

Y_t InterpolationPolynom::operator()(X_t x)
{
    auto l = [this](std::size_t i, const X_t &x)
    {
        X_t li = 1;
        for(std::size_t j = 0; j < m_dots.size(); j++) {
            if (i == j) {
                continue;
            }

            li *= ((x - m_dots[j].x) / (m_dots[i].x - m_dots[j].x));
        }

        return li;
    };

    Y_t y = 0;

    //O = n^2
    for(std::size_t i = 0; i < m_dots.size(); i++) {
        //todo: распаралеливание каждых, например, 5к точек
        y += (l(i, x) * m_dots[i].y);
    }

    return y;
}

CanonicPolynom InterpolationPolynom::toCanonicPolynom() const
{
    int n = m_dots.size();
    if (n <= 0) {
        return CanonicPolynom({});
    }

    // Шаг 1: Вычисление разделённых разностей
    std::vector<std::vector<ValueType>> divDiff(n, std::vector<ValueType>(n));
    for (int i = 0; i < n; ++i) {
        divDiff[i][0] = m_dots[i].y;
    }

    for (int j = 1; j < n; ++j) {
        for (int i = j; i < n; ++i) {
            divDiff[i][j] = (divDiff[i][j-1] - divDiff[i-1][j-1]) / (m_dots[i].x - m_dots[i-j].x);
        }
    }

    // Шаг 2: Получение коэффициентов полинома Ньютона (диагональные элементы)
    coefs_t newtonCoeffs(n);
    for (int i = 0; i < n; ++i) {
        newtonCoeffs[i] = divDiff[i][i];
    }

    // Шаг 3: Преобразование в каноническую форму
    coefs_t canonicalCoeffs(n, 0.0);

    canonicalCoeffs[0] = newtonCoeffs[0];

    coefs_t newtonBasisPoly(n, 0.0);
    newtonBasisPoly[0] = 1.0;

    for (int i = 1; i < n; ++i) {
        coefs_t nextBasisPoly(n, 0.0);

        for (int k = 0; k < i; ++k) {
            nextBasisPoly[k+1] += newtonBasisPoly[k];
        }

        for (int k = 0; k < i; ++k) {
            nextBasisPoly[k] -= newtonBasisPoly[k] * m_dots[i-1].x;
        }
        newtonBasisPoly = nextBasisPoly;

        for (int j = 0; j <= i; ++j) {
            canonicalCoeffs[j] += newtonCoeffs[i] * newtonBasisPoly[j];
        }
    }

    return CanonicPolynom(canonicalCoeffs);
}

PartedCanonicPolynom InterpolationPolynom::toPartedCanonicPolynom() const
{
    std::set<dot_t> sortedDots(m_dots.begin(), m_dots.end());

    return PartedCanonicPolynom(sortedDots);
}

ZeroPolynom::ZeroPolynom(const xs_t &xs)
    : m_roots{xs}
{
    m_coefs = CanonicPolynom::coefsFromRoots(xs);
}

PartedCanonicPolynom ZeroPolynom::toPartedCanonicPolynom() const
{
    std::set<dot_t> sortedDots;
    for(const auto &root : m_roots) {
        sortedDots.insert(dot_t{root, 0});
    }

    return PartedCanonicPolynom(sortedDots, false);
}

Range::Range(X_t left, X_t right)
    : m_left{left}
    , m_right{right}
{
    assert(cmp(right, left) != -1);
    //todo: этот assert для заполнения
    //    assert(cmp(right - left + 1, PartedCanonicPolynom::Partition) == 0);
}

bool Range::inRangeStrict(X_t x) const
{
    return cmp(m_right, x) == 1 && cmp(m_left, x) == -1;
}

bool Range::inRange(X_t x) const
{
    return cmp(m_right, x) >= 0 && cmp(m_left, x) <= 0;
}

Range::pos_t Range::isCrossStrict(const Range &other) const
{
    if (*this == other) {
        return equal;
    }

    if (inRange(other.m_left) && inRange(other.m_right)) {
        return inside;
    }

    if (other.inRange(m_left) && other.inRange(m_right)) {
        return outside;
    }

    if (inRangeStrict(other.m_left) || inRangeStrict(other.m_right) ||
        other.inRangeStrict(m_left) || other.inRangeStrict(m_right)) {
        return crossed;
    }

    if(cmp(m_left, other.m_right) >= 0) {
        return right;
    } else {
        return left;
    }
}

Range Range::crossByStrict(const Range &other) const
{
    if (!(isCrossStrict(other) & (crossed | inside | outside))) {
        return {0, 0};
    }

    if(inRangeStrict(other.m_left)  && other.inRangeStrict(m_right)) {
        return {other.m_left, m_right};
    }

    if(other.inRangeStrict(m_left)  && inRangeStrict(other.m_right)) {
        return {m_left, other.m_right};
    }

    return {0, 0};
}

X_t Range::leftBound() const
{
    return m_left;
}

X_t Range::rightBound() const
{
    return m_right;
}

Range Range::fromUnsorted(X_t a, X_t b)
{
    if (cmp(a, b) == -1) {
        return {a, b};
    }

    return {b, a};
}

bool operator<(const Range &a, const Range &b)
{
    return cmp(a.rightBound(), b.leftBound()) <= 0;
}

bool operator==(const Range &a, const Range &b)
{
    return cmp(a.leftBound(), b.leftBound()) == 0 && cmp(a.rightBound(), b.rightBound()) == 0;
}

bool operator<=(const Range &a, const Range &b)
{
    return a < b || a == b;
}

std::ostream &operator<<(std::ostream &out, const Range &r)
{
    return out << "{" << r.leftBound().get_d() << ", " << r.rightBound().get_d() << "}";
}

PartedCanonicPolynom::PartedCanonicPolynom(const map &map)
    : m_map{map}
{
}

PartedCanonicPolynom::PartedCanonicPolynom(const std::set<dot_t> &sortedDots, bool fromInterpolation)
{
    int partition = fromInterpolation ? PartedCanonicPolynom::Partition : PartedCanonicPolynom::Partition - 1;

    for(auto it = sortedDots.begin(); it != sortedDots.end();) {
        dots_t dots;
        for(int i = 0; i < partition; i++) {
            dots.push_back(*it);

            it++;
            if (it == sortedDots.end()) {
                break;
            }
        }
        if (it != sortedDots.end()) {
            it--;
        }

        X_t rangeEnd = dots.back().x;
        for(std::size_t i = 0; i < partition - dots.size(); i++) {
            rangeEnd++;
        }

        if (fromInterpolation) {
            m_map.insert({dots.front().x, rangeEnd}, InterpolationPolynom(dots).toCanonicPolynom());
        } else {
            xs_t xs;
            for (const auto &[x, y] : dots) {
                xs.insert(x);
            }

            m_map.insert({dots.front().x, rangeEnd}, CanonicPolynom(CanonicPolynom::coefsFromRoots(xs)));
        }
    }
}

Y_t PartedCanonicPolynom::operator()(X_t x)
{
    return m_map[x](x);
}

CustomPolynom PartedCanonicPolynom::operator/(PartedCanonicPolynom &other)
{
    RangeMap<CustomPolynom> result;

    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        CanonicPolynom f1 = it->second;
        CanonicPolynom f2 = itOther->second;
        result.insert(currentRange, f1 / f2);
    });

    return CustomPolynom([result](X_t x) mutable {return result[x](x);});
}

PartedCanonicPolynom PartedCanonicPolynom::operator()(const CanonicPolynom &other) const
{
    assert(other.degree() <= 1);

    auto normInterval = [&other](X_t bound) -> X_t
    {
        if (other.degree() > 0 && other[1] == 0) {
            return bound -= other[0];
        }

        other.degree() >= 0 ? bound -= other[0] : 0.0;
        other.degree() > 0 ? bound /= other[1] : 0.0;

        return bound;
    };

    map result;
    for(auto it = m_map.cbegin(); it != m_map.cend(); it++) {
        auto left = normInterval(it->first.leftBound());
        auto right = normInterval(it->first.rightBound());

        result.insert(Range{left, right}, it->second(other));
    }

    return PartedCanonicPolynom(result);
}

PartedCanonicPolynom PartedCanonicPolynom::operator()(const PartedCanonicPolynom &other) const
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
//        auto otherF = itOther->second;
//        std::vector<Range> ranges;

//        switch(otherF.degree()) {
//        case 0: {
//            X_t left = currentRange.leftBound(), right = currentRange.rightBound();

//            for(auto bound : {&left, &right}) {
//                *bound -= otherF[0];
//            }
//            ranges.push_back({left, right});
//            break;
//        }
//        case 1: {
//            X_t left = currentRange.leftBound(), right = currentRange.rightBound();

//            for(auto bound : {&left, &right}) {
//                otherF.degree() >= 0 ? *bound -= otherF[0] : 0.0;
//                otherF.degree() > 0 ? *bound /= otherF[1] : 0.0;
//            }

//            ranges.push_back({left, right});
//            break;
//        }
//        //todo: тестирование
//        case 2: {
//            X_t left = currentRange.leftBound(), right = currentRange.rightBound();

//            std::vector<X_t> roots;
//            for(auto bound : {left, right}) {
//                X_t a = otherF[2], b = otherF[1], c = otherF[0] - bound;

//                auto discriminant = b * b - 4 * a * c;

//                if (cmp(discriminant, 0) >= 0) {
//                    roots.push_back((-b + sqrt(discriminant)) / (2 * a));
//                    roots.push_back((-b - sqrt(discriminant)) / (2 * a));
//                } else {
//                    assert(false);
//                }
//            }

//            auto rangesCount = roots.size() / 2;
//            assert(roots.size() % 2 == 0);

//            for(std::size_t i = 0; i < rangesCount; i++) {
//                auto r = Range::fromUnsorted(roots[i], roots[i + rangesCount]);
//                std::cout << r << std::endl;
//                ranges.push_back(r);
//            }

//            break;
//        }
//        default:
//            assert(false);
//        }

//        for(const auto &r : ranges) {
//            result.insert(r, it->second(otherF));
//        }


        result.insert(currentRange, it->second(itOther->second));
    });

    return PartedCanonicPolynom(result);
}

PartedCanonicPolynom PartedCanonicPolynom::operator+(const PartedCanonicPolynom &other) const
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        result.insert(currentRange, it->second + itOther->second);
    });

    return PartedCanonicPolynom(result);
}

PartedCanonicPolynom PartedCanonicPolynom::operator*(const PartedCanonicPolynom &other) const
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        result.insert(currentRange, it->second * itOther->second);
    });

    return PartedCanonicPolynom(result);
}

PartedCanonicPolynom PartedCanonicPolynom::operator-(const PartedCanonicPolynom &other) const
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        result.insert(currentRange, it->second - itOther->second);
    });

    return PartedCanonicPolynom(result);
}

void PartedCanonicPolynom::operator+=(const PartedCanonicPolynom &other)
{
    if (m_map.size() == 0) {
        m_map = other.m_map;
        return;
    }

    *this = *this + other;
}

void PartedCanonicPolynom::operatorPrivate(const PartedCanonicPolynom &other, operatorPred_t pred) const
{
    auto mapCopy = m_map;
    auto mapCopyOther = other.m_map;

    {
        auto copyFirst = mapCopy.cbegin();
        auto copyOtherFirst = mapCopyOther.cbegin();

        if (cmp(copyFirst->first.leftBound(), copyOtherFirst->first.leftBound()) == -1) {
            mapCopyOther.insert(Range{copyFirst->first.leftBound(), copyOtherFirst->first.leftBound()}, copyOtherFirst->second);
        } else if (cmp(copyFirst->first.leftBound(), copyOtherFirst->first.leftBound()) == 1){
            mapCopy.insert(Range{copyOtherFirst->first.leftBound(), copyFirst->first.leftBound()}, copyFirst->second);
        }
    }
    {
        auto copyLast = (--mapCopy.cend());
        auto copyOtherLast = (--mapCopyOther.cend());

        if (cmp(copyLast->first.rightBound(), copyOtherLast->first.rightBound()) == -1) {
            mapCopy.insert(Range{copyLast->first.rightBound(), copyOtherLast->first.rightBound()}, copyLast->second);
        } else if (cmp(copyLast->first.rightBound(), copyOtherLast->first.rightBound()) == 1) {
            mapCopyOther.insert(Range{copyOtherLast->first.rightBound(), copyLast->first.rightBound()}, copyOtherLast->second);
        }
    }

    Range currentRange = {0, 0};
    auto it = mapCopy.cbegin();
    auto itOther = mapCopyOther.cbegin();

    auto end = mapCopy.cend();
    auto otherEnd = mapCopyOther.cend();

    auto safeIt = [](map::const_iterator it, map::const_iterator end){return it == end ? (--end) : it;};
    auto safeInc = [](map::const_iterator &it, map::const_iterator end){it != end ? it++ : it;};

    auto itRange = [&it, &mapCopy, &safeIt](){return safeIt(it, mapCopy.cend())->first;};
    auto otherRange = [&itOther, &mapCopyOther, &safeIt](){ return safeIt(itOther, mapCopyOther.cend())->first;};

    auto inc = [&it, &mapCopy, &safeInc](){safeInc(it, mapCopy.cend());};
    auto incOther = [&itOther, &mapCopyOther, &safeInc](){safeInc(itOther, mapCopyOther.cend());};

    auto currIt = [&it, &end, &safeIt](){return safeIt(it, end);};
    auto currItOther = [&itOther, &otherEnd, &safeIt](){return safeIt(itOther, otherEnd);};

    auto predCall = [&currIt, &currItOther, &currentRange, &pred](){pred(currIt(), currItOther(), currentRange);};

    auto moveDuringCross = [&safeIt, &safeInc](map::const_iterator *it, map::const_iterator end, Range cross)
    {
        if (cmp(safeIt((*it), end)->first.rightBound(), cross.leftBound()) == 0 ||
            cmp(safeIt((*it), end)->first.leftBound(), cross.rightBound()) == 0 ||
            cmp(safeIt((*it), end)->first.rightBound(), cross.rightBound()) == 0) {
            safeInc(*it, end);
        }
    };

    auto onInside = [&currentRange, &predCall](Range r, Range cross)
    {
        if (cmp(r.leftBound(), cross.leftBound()) == -1 &&
            cmp(currentRange.rightBound(), r.leftBound()) <= 0) {
            currentRange = Range{r.leftBound(), cross.leftBound()};
            predCall();
        }

        currentRange = cross;
        predCall();
    };

    while(it != end || itOther != otherEnd) {
//        std::cout << itRange() << " " << otherRange() << std::endl;
//        if (currentRange.inRange(1318)) {
//            int i = 0;
//            i += 1;;
//        }

        switch(itRange().isCrossStrict(otherRange())) {
        case Range::equal: {
            currentRange = itRange();
            predCall();

            inc();
            incOther();
            break;
        }
        case Range::crossed: {
            auto cross = itRange().crossByStrict(otherRange());

            map::const_iterator *left, *right, leftEnd, rightEnd;
            if (cmp(itRange().leftBound(),  otherRange().leftBound()) == -1) {
                left = &it;
                leftEnd = end;
                right = &itOther;
                rightEnd = otherEnd;
            } else {
                left = &itOther;
                leftEnd = otherEnd;
                right = &it;
                rightEnd = end;
            }

            if (cmp(currentRange.rightBound(), (*left)->first.leftBound()) <= 0) {
                currentRange = Range{(*left)->first.leftBound(), cross.leftBound()};
                predCall();
            }

            moveDuringCross(left, leftEnd, cross);
            moveDuringCross(right, rightEnd, cross);

            currentRange = cross;
            predCall();

            break;
        }
        case Range::inside: {
            auto cross = otherRange();

            onInside(itRange(), cross);

            moveDuringCross(&it, end, cross);
            moveDuringCross(&itOther, otherEnd, cross);

            break;
        }
        case Range::outside: {
            auto cross = itRange();

            onInside(otherRange(), cross);

            moveDuringCross(&it, end, cross);
            moveDuringCross(&itOther, otherEnd, cross);

            break;
        }

        case Range::left: {
            if (it == end) {
                goto right;
            }

        left:
            currentRange = itRange();
            predCall();

            inc();
            break;
        }
        case Range::right: {
            if (itOther == otherEnd) {
                goto left;
            }

        right:
            currentRange = otherRange();
            predCall();

            incOther();
            break;
        }
        }
    }
}

const int PartedCanonicPolynom::Partition = 3;

}
