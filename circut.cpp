#include "circut.h"

namespace snrk {

template<typename V>
Value<V>::Value(V v)
    : m_data{std::make_shared<V>(new V(v))}
{

}

template<typename V>
bool Value<V>::operator==(const Value &other)
{
    return m_data == other.m_data;
}

template<typename V>
Gate<V>::Gate(type_t type, input_t input, value_t output)
    : m_type{type}
    , m_input{input}
    , m_output{output}
{

}

template<typename V>
Circut<V>::Circut(const input_t &inputX, const input_t &inputW)
    : m_inputX{inputX}
    , m_inputW{inputW}
{

}

template<typename V>
std::size_t Circut<V>::size() const
{
    return m_gates.size();
}

template<typename V>
std::size_t Circut<V>::inputSize() const
{
    return m_inputX.size() + m_inputW.size();
}

template<typename V>
std::size_t Circut<V>::degree() const
{
    return 3 * size() + inputSize();
}

template<typename V>
void Circut<V>::addGate(const gate_t &gate)
{
    /*todo: мютекс, если многопоточка*/

    m_gates.push_back(gate);
}

}
