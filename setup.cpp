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
            t.m_map[*cw++] = elment;
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

S S::generate(const Circut &circut)
{
    S s;

    for(std::size_t i = 0; i < circut.size(); i++) {
        s.m_map[i] = circut.m_gates[i].m_type;
    }

    return s;
}

Gate::type_t S::operator()(std::size_t gateNum)
{
    return m_map.at(gateNum);
}

W W::generate(const witnesses_t &witnesses, const Circut &circut)
{
    W w;
    std::map<value_t, std::shared_ptr<cond_t>> duplicates;
    auto insert = [&duplicates](value_t key, witness_t value)
    {
        if (!duplicates[key]) {
            duplicates[key] = std::make_shared<cond_t>();
        }

        duplicates[key]->insert(value);
    };

    auto cw = witnesses.cbegin();

    auto fillMap = [&cw, &duplicates, &insert](const values_t &row)
    {
        for(const auto &elment : row) {
            insert(elment, *cw++);
        }
    };

    fillMap(circut.m_inputX);
    fillMap(circut.m_inputW);

    for(const auto& gate : circut.m_gates) {
        insert(gate.m_input.a, *cw++);
        insert(gate.m_input.b, *cw++);
        insert(gate.m_output, *cw++);
    }

    for(const auto &[_, condition] : duplicates) {
        for(const auto &witness : *condition) {
            w.m_map[witness] = condition;
        }
    }

    return w;
}

W::cond_t W::operator()(witness_t w) const
{
    return *m_map.at(w);
}

GlobalParams setup(const Circut &circut)
{
    witnesses_t witnesses;
    witnesses.resize(circut.degree());

    witness_t wGenerator = -circut.inputSize();
    std::generate(witnesses.begin(), witnesses.end(),
    [&wGenerator]() -> witness_t
    {
        return wGenerator++;
    });

    /*кадый в отдельный поток, т.к. читаем witness и circut*/
    T t = T::generate(witnesses, circut);
    S s = S::generate(circut);
    W w = W::generate(witnesses, circut);

    return {.pp = {t, s, w}, .vp = {0, 0, 0}};
}

}
