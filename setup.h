#ifndef SETUP_H
#define SETUP_H

#include "circut.h"
#include <unordered_map>
#include <utility>
#include <unordered_set>

namespace snrk {

typedef int witness_t;

template <typename V>
class T {
public:
    static T generate(const std::list<witness_t> &witnesses, const Circut<V> &circut);

    V operator()(witness_t w) const;

private:
    T() = default;

    std::unordered_map<witness_t, V> m_map;
};

template <typename V>
class W {
    using cond_t = std::unordered_set<witness_t>;
public:
    static W generate(const std::list<witness_t> &witnesses, const Circut<V> &circut);

    /*todo: подумать, как проверить иначе T(y) = T(W(y))*/
    cond_t operator()(witness_t w) const;

private:
    std::list<cond_t> m_conditions;

};

template <typename V>
struct GlobalParams {
    using ProverParams = std::pair<T<V>, W<V>>;
    using VerifierParams = std::pair<int, int>; /*todo: func commit*/

    ProverParams pp;
    VerifierParams vp;
};

template <typename V>
GlobalParams<V> setup(const Circut<V> &circut);

}

#endif // SETUP_H
