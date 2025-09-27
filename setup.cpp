#include "setup.h"

#include <algorithm>
#include <map>

namespace snrk {

T T::generate(const witnesses_t &witnesses, const Circut &circut)
{
    T t;

    auto cw = witnesses.cbegin();
    auto fillMap = [&cw, &t](const values_t &row)
    {
        for(const auto &elment : row) {
            t.m_map[*cw] = elment;
            cw++;
        }
    };

    fillMap(circut.m_inputX);
    fillMap(circut.m_inputW);

    for(const auto& gate : circut.m_gates) {
        t.m_map[*cw++] = gate.m_input.a;
        t.m_map[*cw++] = gate.m_input.b;
        t.m_map[*cw++] = gate.m_output;
    }

    return t;
}

value_t T::operator()(witness_t w) const
{
    return m_map.at(w);
}

W W::generate(const witnesses_t &witnesses, const Circut &circut)
{
    W w;
    std::map<value_t, cond_t> duplicates;

    auto cw = witnesses.cbegin();

    auto fillMap = [&cw, &duplicates](const values_t &row)
    {
        for(const auto &elment : row) {
            duplicates[elment].insert(*cw);
            cw++;
        }
    };

    fillMap(circut.m_inputX);
    fillMap(circut.m_inputW);

    for(const auto& gate : circut.m_gates) {
        duplicates[gate.m_input.a].insert(*cw++);
        duplicates[gate.m_input.b].insert(*cw++);
        duplicates[gate.m_output].insert(*cw++);
    }

    for(const auto &[_, condition] : duplicates) {
        w.m_conditions.push_back(condition);
    }

    return w;
}

W::cond_t W::operator()(witness_t w) const
{
    /*todo: тут линейный поиск*/
    auto founded = std::find_if(m_conditions.begin(), m_conditions.end(),
    [&w](const cond_t &condition)
    {
        return condition.count(w) > 0;
    });

    if (founded == m_conditions.end())
    {
        return {};
    }

    return *founded;
}

GlobalParams setup(const Circut &circut)
{
    std::list<witness_t> witnesses;
    witnesses.resize(circut.degree());

    witness_t wGenerator = -circut.inputSize();
    std::generate(witnesses.begin(), witnesses.end(),
    [&wGenerator]() -> witness_t
    {
        return wGenerator++;
    });

    T t = T::generate(witnesses, circut);
    W w = W::generate(witnesses, circut);

    return {.pp = {t, w}, .vp = {0, 0}};
}

}
