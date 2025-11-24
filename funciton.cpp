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

    // Проверка на деление на ноль-полином или степень 0
//    if (nOther < 0 || (nOther == 0 && std::abs(other.m_coefs[0]) < 1e-9)) {
//        // Обработка ошибки (например, выброс исключения)
//        throw std::invalid_argument("Division by zero polynomial or zero constant.");
//    }

    // Если степень делимого меньше степени делителя, частное равно 0, остаток равен делимому
    if (n < nOther) {
        res = CanonicPolynom(1); // Полином 0 степени (константа 0.0)
        res[0] = 0.0;
        rem = *this; // Остаток - исходный полином
        return CustomPolynom::generate([rem, other](X_t x) mutable -> Y_t
        {
            return (rem(x) == 0 ? Y_t(0) : rem(x) / other(x));
        });
    }

    // Инициализируем остаток как копию текущего полинома (делимого)
    // Мы будем модифицировать этот остаток в процессе деления
    coefs_t remainder_coefs = m_coefs;

    // Определяем степень итогового частного и инициализируем его
    int degree_res = n - nOther;
    res = CanonicPolynom(degree_res + 1); // Устанавливаем правильный размер
    for(int i = 0; i <= degree_res; ++i) {
        res[i] = 0.0;
    }

    // Основной цикл длинного деления
    for (int i = degree_res; i >= 0; --i) {
        // Коэффициент текущего члена частного
        ValueType factor = remainder_coefs[i + nOther] / other.m_coefs[nOther];
        res[i] = factor;

        // Вычитаем (factor * x^i * other) из текущего остатка
        for (int j = 0; j <= nOther; ++j) {
            remainder_coefs[i + j] -= factor * other.m_coefs[j];
        }
    }

    // remainder_coefs теперь содержит коэффициенты остатка.
    // Старшие члены (с i + nOther >= nOther) должны быть близки к нулю.
    // Нам нужно скопировать только значимые члены остатка в rem.

    // Остаток имеет степень не более nOther - 1
    rem = CanonicPolynom(nOther); // Размер nOther для индексов 0 .. nOther-1

    for(int j = 0; j < nOther; ++j) {
         // Используем оставшиеся коэффициенты
         rem[j] = remainder_coefs[j];
    }

    return CustomPolynom::generate([res, rem, other](X_t x) mutable -> Y_t
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

    return PartedCanonicPolynom::generate(sortedDots);
}

ZeroPolynom ZeroPolynom::generate(const xs_t &xs)
{
    ZeroPolynom z;

    z.m_coefs = CanonicPolynom::coefsFromRoots(xs);
    z.m_roots = xs;

    return z;
}

PartedCanonicPolynom ZeroPolynom::toPartedCanonicPolynom() const
{
    std::set<dot_t> sortedDots;
    for(const auto &root : m_roots) {
        sortedDots.insert(dot_t{root, 0});
    }

    return PartedCanonicPolynom::generate(sortedDots, false);
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

PartedCanonicPolynom PartedCanonicPolynom::generate(map map)
{
    PartedCanonicPolynom p;

    p.m_map = map;

    return p;
}

PartedCanonicPolynom PartedCanonicPolynom::generate(std::set<dot_t> sortedDots, bool fromInterpolation)
{
    int partition = fromInterpolation ? PartedCanonicPolynom::Partition : PartedCanonicPolynom::Partition - 1;

    PartedCanonicPolynom::map map;
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
            map.insert({dots.front().x, rangeEnd}, InterpolationPolynom::generate(dots).toCanonicPolynom());
        } else {
            xs_t xs;
            for (const auto &[x, y] : dots) {
                xs.insert(x);
            }

            map.insert({dots.front().x, rangeEnd}, CanonicPolynom::generate(CanonicPolynom::coefsFromRoots(xs)));
        }
    }

    return generate(map);
}

Y_t PartedCanonicPolynom::operator()(X_t x)
{
    return m_map[x](x);
}

CustomPolynom PartedCanonicPolynom::operator/(PartedCanonicPolynom &other)
{
    RangeMap<CustomPolynom> result;
    //todo: assert(m_map == other.m_map);

    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        CanonicPolynom f1 = it->second;
        CanonicPolynom f2 = itOther->second;
        std::cout << "left: " << it->first << " | " << "right: " << itOther->first << std::endl;
        ValueType mid = (currentRange.leftBound() + currentRange.rightBound()) / 2;
        std::cout << currentRange << " " << f1(currentRange.leftBound()) << " ? " << f1(currentRange.leftBound()) << std::endl;
        result.insert(currentRange, f1 / f2);
    });

    return CustomPolynom::generate([result](X_t x) mutable {return result[x](x);});
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

    return PartedCanonicPolynom::generate(result);
}

PartedCanonicPolynom PartedCanonicPolynom::operator()(const PartedCanonicPolynom &other) const
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
        auto otherF = itOther->second;
        assert(otherF.degree() <= 1); //todo: решение уравнения 2й степени

        auto normInterval = [&other = otherF](X_t bound) -> X_t
        {
            if (other.degree() > 0 && other[1] == 0) {
                return bound -= other[0];
            }

            other.degree() >= 0 ? bound -= other[0] : 0.0;
            other.degree() > 0 ? bound /= other[1] : 0.0;

            return bound;
        };

        auto left = normInterval(currentRange.leftBound());
        auto right = normInterval(currentRange.rightBound());

        result.insert(Range{left, right}, it->second(otherF));
    });

    return PartedCanonicPolynom::generate(result);
}

PartedCanonicPolynom PartedCanonicPolynom::operator+(const PartedCanonicPolynom &other) const
{
    map result;
    operatorPrivate(other, [&result](map::const_iterator it, map::const_iterator itOther, Range currentRange)
    {
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
