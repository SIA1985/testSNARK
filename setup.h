#ifndef SETUP_H
#define SETUP_H

#include "circut.h"
#include "funciton.h"

#include <unordered_map>
#include <utility>
#include <unordered_set>
#include <algorithm>

namespace snrk { 

typedef int witness_t;
typedef std::vector<witness_t> witnesses_t;

typedef Lagrange T_t;
typedef Lagrange S_t;

class W_t {
    using cond_t = std::unordered_set<witness_t>;
public:
    static W_t generate(const witnesses_t &witnesses, const Circut &circut);

    cond_t operator()(witness_t w) const;

private:
    std::unordered_map<witness_t, std::shared_ptr<cond_t>> m_map;

};

class GlobalParams {
    using ProverParams_t = struct{T_t t; S_t s; W_t w;};
    using VerifierParams_t = struct{Y_t comT; Y_t comS; int comW;};
    using TG_t = struct{value_t t; int G;};

public:
    GlobalParams(const Circut &circut);

    TG_t TG();

    ProverParams_t PP();

    VerifierParams_t VP();

private:
    void generateT(const witnesses_t &witnesses, const Circut &circut);

    void generateS(const Circut &circut);

    void generateW(const witnesses_t &witnesses, const Circut &circut);

    TG_t m_TG;

    T_t m_T;
    S_t m_S;
    W_t m_W;
};

}

#endif // SETUP_H
