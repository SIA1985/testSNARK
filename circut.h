#ifndef CIRCUT_H
#define CIRCUT_H

#include <list>

template <typename V>
struct Circut
{
    typedef std::list<V>    row;
    typedef std::list<row>  table;

public:
    std::size_t size() const;
    std::size_t inputSize() const;
    std::size_t degree() const;

    row m_input_x{};
    row m_input_w{};

    table m_gates{};
};

#endif // CIRCUT_H
