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

class PolynomSubstitutionProof : public Proof
{
    using commit_t = X_t;

public:
    PolynomSubstitutionProof(commit_t comF, commit_t comQ, dot_t toProve);

    void setGp(X_t t, int G);

    virtual bool check() override;

private:
    commit_t m_comF;
    commit_t m_comQ;

    dot_t m_toProve;

    X_t m_t;
    int m_G;
};

}

#endif // PROOF_H
