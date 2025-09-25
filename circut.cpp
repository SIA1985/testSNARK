#include "circut.h"

snrk::Gate::Gate(type_t type, input_t input, value_t output)
    : m_type{type}
    , m_input{input}
    , m_output{output}
{

}

snrk::Circut::Circut(const input_t &inputX, const input_t &inputW)
    : m_inputX{inputX}
    , m_inputW{inputW}
{

}

std::size_t snrk::Circut::size() const
{
    return m_gates.size();
}

std::size_t snrk::Circut::inputSize() const
{
    return m_inputX.size() + m_inputW.size();
}

std::size_t snrk::Circut::degree() const
{
    return 3 * size() + inputSize();
}

void snrk::Circut::addGate(const gate_t &gate)
{
    /*todo: мютекс, если многопоточка*/

    m_gates.push_back(gate);
}
