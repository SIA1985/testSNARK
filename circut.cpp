#include "circut.h"

#include <cmath>


namespace snrk {

Value::Value()
    : m_value{std::make_shared<ValueType>()}
{
}

Value::Value(ValueType value)
    : m_value{std::make_shared<ValueType>(value)}
{
}

Value::operator int() const
{
    return static_cast<int>(m_value->get_si());
}

Value::operator double() const
{
    return m_value->get_d();
}

Value::operator GateType_t() const
{
    if (floor(*m_value) != *m_value) {
        return Unknown;
    }

    return static_cast<GateType_t>(m_value->get_ui());
}

Value::operator ValueType() const
{
    return *m_value;
}

ValueType* Value::get() const
{
    return m_value.get();
}

bool operator<(const Value &a, const Value &b)
{
    return *a.m_value < *b.m_value;
}

//bool operator==(const Value &a, const Value &b)
//{
//    return a.m_value.get() == b.m_value.get();
//}


Gate::Gate(GateType_t type, input_t input, value_t output)
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
