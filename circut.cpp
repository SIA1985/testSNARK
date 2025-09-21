#include "circut.h"


template<typename V>
std::size_t Circut<V>::size() const
{
    return m_gates.size();
}

template<typename V>
std::size_t Circut<V>::inputSize() const
{
    return m_input_x.size() + m_input_w.size();
}

template<typename V>
std::size_t Circut<V>::degree() const
{
    return 3 * size() + inputSize();
}
