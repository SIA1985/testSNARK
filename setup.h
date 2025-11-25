#ifndef SETUP_H
#define SETUP_H

#include "circut.h"
#include "funciton.h"

#include <unordered_map>
#include <utility>
#include <unordered_set>
#include <algorithm>

namespace snrk { 

typedef double witness_t;
typedef std::vector<witness_t> witnesses_t;

witnesses_t genWitnesses(witness_t start, std::size_t count);

typedef InterpolationPolynom T_t;
typedef InterpolationPolynom S_t;

class W_t {
    using cond_t = std::unordered_set<witness_t>;
public:
    static W_t generate(const witnesses_t &witnesses, const Circut &circut);

    cond_t operator()(witness_t w) const;

private:
    std::unordered_map<witness_t, std::shared_ptr<cond_t>> m_map;
};

struct SplittedT_t
{
    InterpolationPolynom left;
    InterpolationPolynom right;
    InterpolationPolynom result;
};

class GlobalParams {
    using ProverParams_t = struct{T_t t; SplittedT_t splittedT; S_t s; W_t w;};

public:
    GlobalParams(const Circut &circut);

    witnesses_t witnesses() const;

    TG_t TG();

    ProverParams_t PP();

private:
    void generateT(const Circut &circut);

    void generateS(const Circut &circut);

    void generateW(const Circut &circut);

    witnesses_t m_witnesses;

    TG_t m_TG;

    T_t m_T;
    SplittedT_t m_splittedT;
    S_t m_S;
    W_t m_W;
};

}

#endif // SETUP_H
