#ifndef CIRCUT_H
#define CIRCUT_H

#include <list>
#include <memory>

template <typename V>
struct Value {
    using D = std::shared_ptr<V>;


    Value(V);

    bool operator==(const Value &other);

private:
    D m_data;
};

template <typename V>
class Gate {
    enum type_t : char {
        Sum = 0,
        Product,
    };

    using value_t = Value<V>;
    using input_t = std::list<value_t>;

public:
    Gate(type_t type, input_t input, value_t output);

private:
    type_t m_type{};
    input_t m_input{};
    value_t m_output{};
};

template <typename V>
struct Circut
{
    using input_t = std::list<V>;
    using gates_t = std::list<Gate<V>>;

public:
    std::size_t size() const;
    std::size_t inputSize() const;
    std::size_t degree() const;

    input_t m_inputX{};
    input_t m_inputW{};

    gates_t m_gates{};
};

#endif // CIRCUT_H
