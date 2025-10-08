#ifndef CIRCUT_H
#define CIRCUT_H

#include <vector>
#include <memory>
#include <any>

namespace snrk {

#define ValueType double

class Value
{
public:
    Value();

    Value(ValueType value);

    operator int() const;
    operator double() const;

    bool operator<(const Value &other);
    bool operator==(const Value &other);

private:
    std::shared_ptr<ValueType> m_value{};
};

bool operator<(const Value &a, const Value &b);
bool operator==(const Value &a, const Value &b);

typedef Value value_t;
typedef std::vector<value_t> values_t;

class Gate {
    using input_t = struct{value_t a; value_t b;};

public:
    enum type_t : char {
        Sum = 0,
        Product,
    };

    Gate(type_t type, input_t input, value_t output);

private:
    type_t m_type{};
    input_t m_input{};
    value_t m_output{0};

    friend class GlobalParams;
    friend class W_t;
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
    friend class W_t;
};

}

#endif // CIRCUT_H
