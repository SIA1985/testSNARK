#ifndef PROOF_H
#define PROOF_H

#include "funciton.h"

namespace snrk {

class Proof
{
public:
    Proof() = default;
    virtual ~Proof() = default;

    virtual bool check() = 0;
};

template <typename X, typename Y>
class PolynomSubstitutionProof : public Proof
{
    using commit_t = X;
    using toProve_t = struct{X u; Y v;};

public:
    PolynomSubstitutionProof(commit_t comF, commit_t comQ, toProve_t toProve)
        : m_comF{comF}
        , m_comQ{comQ}
        , m_toProve{toProve}
    {
    }

    void setGp(X t, int G)
    {
        m_t = t;
        m_G = G;
    }

    virtual bool check() override
    {
        /*todo: сравнение float*/
        return (m_t - m_toProve.u) * m_comQ == m_comF - m_toProve.v * m_G;
    }

private:
    commit_t m_comF;
    commit_t m_comQ;

    toProve_t m_toProve;

    X m_t;
    int m_G;
};

}

#endif // PROOF_H
