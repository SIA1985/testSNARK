#ifndef CIRCUT_H
#define CIRCUT_H

#include <list>
#include <memory>

namespace snrk {

template <typename V>
class Value {
    using D = std::shared_ptr<V>;

public:
    Value(V v);

    bool operator==(const Value &other);

private:
    D m_data;
};

template <typename V>
class Gate {
    using value_t = Value<V>;
    using input_t = std::list<value_t>;

public:
    enum type_t : char {
        Sum = 0,
        Product,
    };

    Gate(type_t type, input_t input, value_t output);

private:
    type_t m_type{};
    input_t m_input{};
    value_t m_output{};
};

template <typename V>
class Circut
{
    using input_t = std::list<V>;
    using gate_t = Gate<V>;
    using gates_t = std::list<Gate<V>>;

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
};

}

#endif // CIRCUT_H
