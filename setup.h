#ifndef SETUP_H
#define SETUP_H

#include <circut.h>
#include <unordered_map>

typedef int W;

template <typename V>
class T {
public:
    static T generate(const Circut<V> &circut);

    V operator()(W w) const;

private:
    T(const Circut<V> &circut);

    std::unordered_map<W, V> m_map;
};

#endif // SETUP_H
