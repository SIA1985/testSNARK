#ifndef SETUP_H
#define SETUP_H

#include "circut.h"
#include <unordered_map>
#include <utility>
#include <unordered_set>
#include <algorithm>

namespace snrk {

typedef int witness_t;
typedef std::list<witness_t> witnesses_t;

class T {
public:
    static T generate(const witnesses_t &witnesses, const Circut &circut);

    value_t operator()(witness_t w) const;

private:
    T() = default;

    std::unordered_map<witness_t, value_t> m_map;
};

class W {
    using cond_t = std::unordered_set<witness_t>;
public:
    static W generate(const witnesses_t &witnesses, const Circut &circut);

    cond_t operator()(witness_t w) const;

private:
    W() = default;

    std::list<cond_t> m_conditions;

};

struct GlobalParams {
    using ProverParams = struct{T t; W w;};
    using VerifierParams = struct{int comT; int comW;}; /*todo: func commit*/

    ProverParams pp;
    VerifierParams vp;
};

GlobalParams setup(const Circut &circut);

}

#endif // SETUP_H
