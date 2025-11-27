#ifndef SETUP_H
#define SETUP_H

#include "circut.h"
#include "polynom.h"

#include <unordered_map>
#include <utility>
#include <unordered_set>
#include <algorithm>

namespace snrk { 

static witness_t wStart = 1.;
witnesses_t genWitnesses(witness_t start, std::size_t count, bool Chebishev = false);

class W_t {
    using cond_t = std::unordered_set<witness_t>;
public:
    W_t() = default;

    W_t(const witnesses_t &witnesses, const Circut &circut);

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
public:
    using opsFromS_t = std::unordered_map<GateType_t, dots_t>;
    using TParams_t = struct{T_t t; SplittedT_t splittedT;};
    using SParams_t = struct{S_t s; opsFromS_t opsFromS;};

    using ProverParams_t = struct{TParams_t TParams; SParams_t SParams; W_t w;};


    GlobalParams(const Circut &circut);

    witnesses_t witnesses() const;
    witnesses_t SWitnesses() const;

    TG_t TG();

    ProverParams_t PP();

private:
    void generateT(const Circut &circut);

    void generateS(const Circut &circut);

    void generateW(const Circut &circut);

    witnesses_t m_witnesses;
    witnesses_t m_SWitnesses;

    TG_t m_TG;

    T_t m_T;
    SplittedT_t m_splittedT;

    S_t m_S;
    opsFromS_t m_opsFromS;

    W_t m_W;
};

}

#endif // SETUP_H
