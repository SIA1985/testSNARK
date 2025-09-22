#ifndef SETUP_H
#define SETUP_H

#include <circut.h>
#include <unordered_map>

typedef int w_t;

template <typename V>
class T {
public:
    static T generate(const Circut<V> &circut);

    V operator()(w_t w) const;

private:
    T() = default;

    std::unordered_map<w_t, V> m_map;
};

template <typename V>
class W {
public:

};

#endif // SETUP_H
