#include "setup.h"

#include <algorithm>

template<typename V>
T<V>::T(const Circut<V> &circut)
{
    using circutRow = typename Circut<V>::row;

    std::list<W> ws;
    ws.resize(circut.degree());

    W w = -circut.inputSize();
    std::generate(ws.begin(), ws.end(),
    [&w]() -> W
    {
        return w++;
    });

    auto cw = ws.cbegin();
    auto fillMap = [&cw, this](const circutRow& row)
    {
        for(const auto &elment : row) {
            m_map[*cw] = elment;
            cw++;
        }
    };

    fillMap(circut.m_input_x);
    fillMap(circut.m_input_w);

    for(const auto& row : circut.m_gates) {
        fillMap(row);
    }
}

template<typename V>
T<V> T<V>::generate(const Circut<V> &circut)
{
    return T(circut);
}

template<typename V>
V T<V>::operator()(W w) const
{
    return m_map.at(w);
}
