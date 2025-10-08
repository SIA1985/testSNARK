#include "circut.h"

namespace snrk {

Value::Value()
    : m_value{std::make_shared<ValueType>()}
{
}

Value::Value(double value)
    : m_value{std::make_shared<ValueType>(value)}
{
}

Value::operator int() const
{
    return static_cast<int>(*m_value);
}

Value::operator double() const
{
    return static_cast<double>(*m_value);
}

bool Value::operator<(const Value &other)
{
    return *m_value < *other.m_value;
}

bool operator<(const Value &a, const Value &b)
{
    return a < b;
}

bool Value::operator==(const Value &other)
{
    return m_value == other.m_value;
}

bool operator==(const Value &a, const Value &b)
{
    return a == b;
}

Gate::Gate(type_t type, input_t input, value_t output)
    : m_type{type}
    , m_input{input}
    , m_output{output}
{

}

Circut::Circut(const input_t &inputX, const input_t &inputW)
    : m_inputX{inputX}
    , m_inputW{inputW}
{

}

std::size_t Circut::size() const
{
    return m_gates.size();
}

std::size_t Circut::inputSize() const
{
    return m_inputX.size() + m_inputW.size();
}

std::size_t Circut::degree() const
{
    return 3 * size() + inputSize();
}

void Circut::addGate(const gate_t &gate)
{
    /*todo: мютекс, если многопоточка*/

    m_gates.push_back(gate);
}

}
