#include "funciton.h"

#include <cmath>
#include <iostream>

namespace snrk {

bool operator<(const dot_t &a, const dot_t &b)
{
    return cmp(a.x, b.x) != 1;
}

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

Y_t CustomPolynom::operator()(X_t x)
{
    return m_customFunction(x);
}

CanonicPolynom CanonicPolynom::generate(coefs_t coefs)
{
    CanonicPolynom c;

    c.m_coefs = coefs;

    return c;
}

CanonicPolynom::coefs_t CanonicPolynom::coefsFromRoots(roots_t roots)
{
    CanonicPolynom result(1); // Начинаем с полинома p(x) = 1
    result[0] = 1.0;

    for (const auto &root : roots) {
        // Создаем полином-двучлен (x - root)
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

    for(j = 0; j <= n; j++) {
            rem[j] = (*this)[j];
            res[j]=0.0;
    }
    for(k = n - nOther; k >= 0; k--) {
        res[k] = rem[nOther + k] / other[nOther];

        for(j = nOther + k - 1; j >= k; j--) {
            rem[j] -= res[k] * other[j - k];
        }
    }
    for(j = nOther; j <= n; j++) {
        rem[j] = 0.0;
    }

    return CustomPolynom::generate([res, rem](X_t x) mutable -> Y_t
    {
        return res(x) + rem(x);
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

    return CanonicPolynom::generate(diff);
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

    return CanonicPolynom::generate(diff);
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

    return CanonicPolynom::generate(result);
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
        return CanonicPolynom::generate({});
    }

    // Начинаем с последнего коэффициента (a_n)
    CanonicPolynom result = CanonicPolynom::generate({m_coefs.back()});

    // Идем с конца до второго элемента (индекс 1, пропуская a_0)
    for (int i = m_coefs.size() - 2; i >= 0; --i) {
        // Умножаем текущий результат на Q(x)
        result *= other;
        // Добавляем следующий коэффициент (a_{n-1}, a_{n-2}, ...)
        result += CanonicPolynom::generate({m_coefs[i]});
    }
    // В конце получаем полный P(Q(x))

    /*
    CanonicPolynom y = CanonicPolynom::generate({m_coefs[0]});
    auto xPow = other;

    for(std::size_t i = 1; i < m_coefs.size(); i++) {
        y += xPow * m_coefs[i];
        xPow *= other;
    }*/


    return result;
}

ValueType &CanonicPolynom::operator[](std::size_t i)
{
    return m_coefs[i];
}

std::size_t CanonicPolynom::degree() const
{
    return m_coefs.size() - 1;
}

CanonicPolynom::CanonicPolynom(std::size_t n)
{
    m_coefs.resize(n, ValueType(0));
}

InterpolationPolynom InterpolationPolynom::generate(const dots_t &dots)
{
    InterpolationPolynom l;

    l.m_dots = dots;

    return l;
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

    for(std::size_t i = 0; i < m_dots.size(); i++) {
        y += (l(i, x) * m_dots[i].y);
    }

    return y;
}

CanonicPolynom InterpolationPolynom::toCanonicPolynom() const
{
    int n = m_dots.size();
    if (n <= 0) {
        return CanonicPolynom::generate({});
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
    std::vector<ValueType> newtonCoeffs(n);
    for (int i = 0; i < n; ++i) {
        newtonCoeffs[i] = divDiff[i][i];
    }

    // Шаг 3: Преобразование в каноническую форму
    std::vector<ValueType> canonicalCoeffs(n, 0.0);

    // Инициализация первого коэффициента
    canonicalCoeffs[0] = newtonCoeffs[0];

    // Вспомогательный полином для раскрытия скобок
    std::vector<ValueType> newtonBasisPoly(n, 0.0);
    newtonBasisPoly[0] = 1.0;

    for (int i = 1; i < n; ++i) {
        // Умножение текущего базисного полинома на (x - x_{i-1})
        std::vector<ValueType> nextBasisPoly(n, 0.0);

        // Член x * newton_basis_poly
        for (int k = 0; k < i; ++k) {
            nextBasisPoly[k+1] += newtonBasisPoly[k];
        }
        // Член -x_{i-1} * newton_basis_poly
        for (int k = 0; k < i; ++k) {
            nextBasisPoly[k] -= newtonBasisPoly[k] * m_dots[i-1].x;
        }
        newtonBasisPoly = nextBasisPoly;

        // Добавление к каноническому полиному
        for (int j = 0; j <= i; ++j) {
            canonicalCoeffs[j] += newtonCoeffs[i] * newtonBasisPoly[j];
        }
    }

    return CanonicPolynom::generate(canonicalCoeffs);
}

PartedCanonicPolynom InterpolationPolynom::toPartedCanonicPolynom() const
{
    std::set<dot_t> sortedDots(m_dots.begin(), m_dots.end());

    PartedCanonicPolynom::map map;
    for(auto it = sortedDots.begin(); it != sortedDots.end();) {
        dots_t dots;
        for(int i = 0; i < PartedCanonicPolynom::Partition; i++) {
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
        for(std::size_t i = 0; i < PartedCanonicPolynom::Partition - dots.size(); i++) {
            rangeEnd++;
        }

        map.insert({dots.front().x, rangeEnd}, InterpolationPolynom::generate(dots).toCanonicPolynom());
    }

    return PartedCanonicPolynom::generate(map);
}

ZeroPolynom ZeroPolynom::generate(const xs_t &xs)
{
    ZeroPolynom z;

    z.m_coefs = coefsFromRoots(xs);

    return z;
}

Range::Range(X_t left, X_t right)
    : m_left{left}
    , m_right{right}
{
    assert(cmp(right, left) != -1);
    assert(cmp(right - left + 1, PartedCanonicPolynom::Partition) == 0);
}

int Range::inRange(X_t x) const
{
    if (cmp(m_left, x) == 1) {
        return -1;
    }
    else if (cmp(m_right, x) == -1) {
        return 1;
    }

    return 0;
}

X_t Range::leftBound() const
{
    return m_left;
}

X_t Range::rightBound() const
{
    return m_right;
}

bool operator<(const Range &a, const Range &b)
{
    return cmp(a.rightBound(), b.leftBound()) <= 0;
}

bool operator==(const Range &a, const Range &b)
{
    return cmp(a.leftBound(), b.leftBound()) == 0 && cmp(a.rightBound(), b.rightBound() == 0);
}

bool operator<=(const Range &a, const Range &b)
{
    return a < b || a == b;
}

PartedCanonicPolynom PartedCanonicPolynom::generate(map map)
{
    PartedCanonicPolynom p;

    p.m_map = map;

    return p;
}

Y_t PartedCanonicPolynom::operator()(X_t x)
{
    return m_map[x](x);
}

PartedCanonicPolynom PartedCanonicPolynom::operator()(const CanonicPolynom &other) const
{
    map result;
    for(auto it = m_map.cbegin(); it != m_map.cend(); it++) {
        result.insert(it->first, it->second(other));
    }

    return PartedCanonicPolynom::generate(result);
}

PartedCanonicPolynom PartedCanonicPolynom::operator()(const PartedCanonicPolynom &other) const
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        result.insert(currentRange, it->second(itOther->second));
    });

    return PartedCanonicPolynom::generate(result);
}

PartedCanonicPolynom PartedCanonicPolynom::operator+(const PartedCanonicPolynom &other) const
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        auto f1 = it->second;
        auto f2 = itOther->second;
        std::cout << currentRange.leftBound() << " " << currentRange.rightBound() << " : " << f1(1) << " " << f2(1) << std::endl;
        result.insert(currentRange, it->second + itOther->second);
    });

    return PartedCanonicPolynom::generate(result);
}

PartedCanonicPolynom PartedCanonicPolynom::operator*(const PartedCanonicPolynom &other) const
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        result.insert(currentRange, it->second * itOther->second);
    });

    return PartedCanonicPolynom::generate(result);
}

CustomPolynom PartedCanonicPolynom::operator/(CanonicPolynom &other)
{
    RangeMap<CustomPolynom> result;
    //todo: assert(m_map == other.m_map);

    for(auto it = m_map.cbegin(); it != m_map.cend(); it++) {
        result.insert(it->first, m_map[it->first] / other);
    }

    return CustomPolynom::generate([result](X_t x) mutable {return result[x](x);});
}

PartedCanonicPolynom PartedCanonicPolynom::operator-(const PartedCanonicPolynom &other) const
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        result.insert(currentRange, it->second - itOther->second);
    });

    return PartedCanonicPolynom::generate(result);
}

void PartedCanonicPolynom::operator+=(const PartedCanonicPolynom &other)
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        result[currentRange] += itOther->second;
    });

    m_map = result;
}

void PartedCanonicPolynom::operatorPrivate(const PartedCanonicPolynom &other, operatorPred_t pred) const
{
    PartedCanonicPolynom left, right;
    bool meLeft = m_map.cbegin()->first.leftBound() < other.m_map.cbegin()->first.leftBound();

    if (meLeft) {
        left = *this;
        right = other;
    } else {
        left = other;
        right = *this;
    }

    Range currentRange = left.m_map.cbegin()->first;
    auto itLeft = left.m_map.cbegin();
    auto itRight = right.m_map.cbegin();

    auto leftVal = [&left, &itLeft](){ return itLeft == left.m_map.cend() ? --left.m_map.cend() : itLeft;};
    auto rightVal = [&right, &itRight](){ return itRight == right.m_map.cend() ? --right.m_map.cend() : itRight;};

    auto leftEnd = left.m_map.cend();
    auto rightEnd = right.m_map.cend();

    while(itLeft != left.m_map.cend() || itRight != right.m_map.cend()) {
        if (meLeft) {
            pred(leftVal(), rightVal(), currentRange);
        } else {
            pred(rightVal(), leftVal(), currentRange);
        }

        if (itLeft != leftEnd) {
            itLeft++;
            currentRange = leftVal()->first;
        } else if (itRight != rightEnd) {
            itRight++;
            currentRange = rightVal()->first;
        }

        //todo: отсечь пересечения
//        if (rightVal()->first == leftVal()->first && itRight != rightEnd) {
//            itRight++;
//        }
    }
}

const int PartedCanonicPolynom::Partition = 3;

}
