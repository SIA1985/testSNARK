#ifndef PROOF_H
#define PROOF_H

#include "funciton.h"
#include "setup.h"

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
    using ptr_t = std::shared_ptr<PolynomSubstitutionProof>;

public:
    static ptr_t forProver(Polynom &f, dot_t toProve, TG_t tG);
    static ptr_t forVerifier(commit_t comF, commit_t comQ, dot_t toProve, TG_t tG);

    virtual bool check() override;

private:
    PolynomSubstitutionProof() = default;

    commit_t m_comF;
    commit_t m_comQ;

    dot_t m_toProve;

    TG_t m_tG;
};

class ZeroTestProof : public Proof
{
    using commit_t = X_t;
    using ptr_t = std::shared_ptr<ZeroTestProof>;

public:
    /*todo: откуда брать r?*/
    /*check: f = p -> f - p = 0*/
    static ptr_t forProver(Polynom &g, Polynom &p, TG_t tG);
    static ptr_t forVerifier(commit_t comF, commit_t comQ, Y_t f_r, Y_t q_r);

    virtual bool check() override;

private:
    ZeroTestProof() = default;

    commit_t m_comF;
    commit_t m_comQ;
};

}

#endif // PROOF_H
