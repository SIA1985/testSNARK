#include "setup.h"

#include <algorithm>

template<typename V>
T<V> T<V>::generate(const Circut<V> &circut)
{
    T t;

    using circutRow = typename Circut<V>::row;

    std::list<w_t> ws;
    ws.resize(circut.degree());

    w_t w = -circut.inputSize();
    std::generate(ws.begin(), ws.end(),
    [&w]() -> w_t
    {
        return w++;
    });

    auto cw = ws.cbegin();
    auto fillMap = [&cw, &t](const circutRow& row)
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
V T<V>::operator()(w_t w) const
{
    return m_map.at(w);
}
