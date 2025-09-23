#include "setup.h"

#include <algorithm>

namespace snrk {

template<typename V>
T<V> T<V>::generate(const std::list<witness_t> &witnesses, const Circut<V> &circut)
{
    T t;

    auto cw = witnesses.cbegin();
    auto fillMap = [&cw, &t](const std::list<V> &row)
    {
        for(const auto &elment : row) {
            t.m_map[*cw] = elment;
            cw++;
        }
    };

    fillMap(circut.m_inputX);
    fillMap(circut.m_inputW);

    for(const auto& gate : circut.m_gates) {
        fillMap(gate.m_input);
        t.m_map[*cw] = gate.m_output;
        cw++;
    }

    return t;
}

template<typename V>
V T<V>::operator()(witness_t w) const
{
    return m_map.at(w);
}

template<typename V>
W<V> W<V>::generate(const std::list<witness_t> &witnesses, const Circut<V> &circut)
{

}

template<typename V>
GlobalParams<V> setup(const Circut<V> &circut)
{
    using T_t = T<V>;
    using W_t = W<V>;

    std::list<witness_t> witnesses;
    witnesses.resize(circut.degree());

    witness_t wGenerator = -circut.inputSize();
    std::generate(witnesses.begin(), witnesses.end(),
    [&wGenerator]() -> witness_t
    {
        return wGenerator++;
    });

    T_t t = T_t::generate(witnesses, circut);
    W_t w = W_t::generate(witnesses, circut);

    return {.pp = {t, w}, .vp = {0, 0}};
}

}
