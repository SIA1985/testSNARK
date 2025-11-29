#ifndef CIRCUT_H
#define CIRCUT_H

#include <vector>
#include <memory>

#include "types.h"

namespace snrk {

class Value
{
public:
    Value();

    Value(ValueType value);

    operator int() const;
    operator double() const;
    operator GateType_t() const;
    operator ValueType() const;

private:
    std::shared_ptr<ValueType> m_value{};

    friend bool operator<(const Value &a, const Value &b);
};

bool operator<(const Value &a, const Value &b);
//bool operator==(const Value &a, const Value &b);

class Gate {

    using input_t = struct{value_t a; value_t b;};
public:
    Gate(GateType_t type, input_t input, value_t output);

private:
    GateType_t m_type{};
    input_t m_input{};
    value_t m_output{0};

    friend class GlobalParams;
};

class Circut
{
    using input_t = std::vector<value_t>;
    using gate_t = Gate;
    using gates_t = std::vector<Gate>;

public:
    Circut(const input_t &inputX, const input_t &inputW);

    std::size_t size() const;

    std::size_t inputSize() const;

    std::size_t degree() const;


    void addGate(const gate_t &gate);

private:
    input_t m_inputX{};
    input_t m_inputW{};

    gates_t m_gates{};

    friend class GlobalParams;
};

}

#endif // CIRCUT_H
