#include "circut.h"

namespace snrk {

template<typename V>
Value<V>::Value(V)
    : m_data{std::make_shared<V>(new V)}
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

}
