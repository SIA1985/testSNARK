#include "funciton.h"

#include <cmath>

namespace snrk {

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

Y_t CustomPolynom::operator()(const X_t &x)
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
    CanonicPolynom result(0 + 1); // Начинаем с полинома p(x) = 1
    result[0] = 1.0;

    for (const auto &root : roots) {
        // Создаем полином-двучлен (x - root)
        CanonicPolynom binomial(1 + 1);
        binomial[0] = -root;
        binomial[1] = 1.0;

        // Умножаем текущий результат на новый двучлен
        result = result * binomial;
    }

    return result.m_coefs;
}

Y_t CanonicPolynom::operator()(const X_t &x)
{
    Y_t y = 0;

    for(std::size_t i = 0; i < m_coefs.size(); i++) {
        y += m_coefs[i] * std::pow(x, i);
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
        res[k] = rem[nOther+k] / other[nOther];

        for(j = nOther + k - 1; j >= k; j--) {
            rem[j] -= res[k] * other[j-k];
        }
    }
    for(j = nOther; j <= n; j++) {
        rem[j] = 0.0;
    }

    return CustomPolynom::generate([res, rem](X_t x) mutable -> Y_t
    {
        auto a = res(x);
        auto b = rem(x);
        return a;
    });
}

/*todo: рефакторинг other = -1 * other*/
/*this - other*/
CanonicPolynom CanonicPolynom::operator-(CanonicPolynom &other)
{
    std::size_t n = m_coefs.size();
    std::size_t nOther = other.m_coefs.size();
    coefs_t diff;

    if (n > nOther) {
        diff.resize(n, 0);

        for(std::size_t i = 0; i < nOther; i++) {
            diff[i] = m_coefs[i] - other.m_coefs[i];
        }
        for(std::size_t i = nOther; i < n; i++) {
            diff[i] = m_coefs[i];
        }
    } else {
        diff.resize(nOther, 0);

        for(std::size_t i = 0; i < n; i++) {
            diff[i] = m_coefs[i] - other.m_coefs[i];
        }
        for(std::size_t i = n; i < nOther; i++) {
            diff[i] = -other.m_coefs[i];
        }
    }

    return CanonicPolynom::generate(diff);
}

CanonicPolynom CanonicPolynom::operator*(CanonicPolynom &other)
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
    auto newtonCoefs = [](const dots_t& points) -> coefs_t {
        int n = points.size();
        coefs_t f(n);
        for (int i = 0; i < n; ++i) {
            f[i] = points[i].y;
        }

        coefs_t coefficients(n);
        coefficients[0] = f[0];

        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < n - i; ++j) {
                f[j] = (f[j + 1] - f[j]) / (points[j + i].x - points[j].x);
            }
            coefficients[i] = f[0];
        }
        return coefficients;
    };

    InterpolationPolynom l;

    l.m_dots = dots;
    l.m_newtonCoefs = newtonCoefs(dots);

    return l;
}

Y_t InterpolationPolynom::operator()(const X_t &x)
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

CanonicPolynom InterpolationPolynom::toClassicPolynom() const
{
    int n = m_newtonCoefs.size();
    coefs_t canonicalCoeffs(n, 0.0);
    coefs_t p(n, 0.0);
    coefs_t prev(n, 0.0);

    // Коэффициент при x^0
    p[0] = 1.0;
    canonicalCoeffs[0] = m_newtonCoefs[0] * p[0];

    // Вычисление остальных коэффициентов
    for (int i = 1; i < n; ++i) {
        prev = p;
        p[0] = -m_dots[i - 1].x * prev[0];
        canonicalCoeffs[0] += m_newtonCoefs[i] * p[0];

        for (int j = 1; j <= i; ++j) {
            p[j] = prev[j - 1] - m_dots[ i - 1].x * prev[j];
            canonicalCoeffs[j] += m_newtonCoefs[i] * p[j];
        }
    }

    return CanonicPolynom::generate(canonicalCoeffs);
}

ZeroPolynom ZeroPolynom::generate(const xs_t &xs)
{
    ZeroPolynom z;

    z.m_coefs = coefsFromRoots(xs);

    return z;
}

/*
 * @brief Вычисляет значение полинома Ньютона в заданной точке.
 *
 * @param newtonCoeffs Вектор коэффициентов Ньютона.
 * @param points Вектор узловых точек, использованных для построения полинома.
 * @param x_eval Точка, в которой нужно вычислить значение полинома.
 * @return Значение полинома в точке x_eval.
 *
double evaluateNewtonPolynomial(const std::vector<double>& newtonCoeffs, const std::vector<Point>& points, double x_eval) {
    int n = newtonCoeffs.size();
    double result = newtonCoeffs[n - 1];

    for (int i = n - 2; i >= 0; --i) {
        result = result * (x_eval - points[i].x) + newtonCoeffs[i];
    }

    return result;
}

int main() {
    // Узловые точки
    std::vector<Point> points = {
        {0, 2},
        {5, -2},
        {10, 7}
    };

*/

}
