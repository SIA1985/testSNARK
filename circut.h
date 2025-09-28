#ifndef CIRCUT_H
#define CIRCUT_H

#include <vector>
#include <memory>
#include <any>

namespace snrk {

#define ValueType int

using value_t = std::shared_ptr<ValueType>;
typedef std::vector<value_t> values_t;

#define Value(value) std::make_shared<ValueType>(value);

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

    friend class T;
    friend class W;
    friend class S;
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

    friend class T;
    friend class W;
    friend class S;
};

}

#endif // CIRCUT_H
