#ifndef PROOF_H
#define PROOF_H

#include "funciton.h"
#include "setup.h"

namespace snrk {

typedef std::function<xs_t(std::size_t)> wGenerator_t;
extern wGenerator_t wGeneratorDefault;

class Proof
{
public:
    Proof() = default;
    virtual ~Proof() = default;

    virtual bool check(wGenerator_t wGenerator = wGeneratorDefault) = 0;
};

class PolynomSubstitutionProof : public Proof
{
    using commit_t = X_t;
    using ptr_t = std::shared_ptr<PolynomSubstitutionProof>;

public:
    static ptr_t forProver(Polynom &f, dot_t toProve, TG_t tG);
    static ptr_t forVerifier(commit_t comF, commit_t comQ, dot_t toProve, TG_t tG);

    virtual bool check(wGenerator_t wGenerator = wGeneratorDefault) override;

private:
    PolynomSubstitutionProof() = default;

    commit_t m_comF;
    commit_t m_comQ;

    dot_t m_toProve;

    TG_t m_tG;

    friend class ZeroTestProof;
};

class ZeroTestProof : public Proof
{
    using commit_t = X_t;
    using ptr_t = std::shared_ptr<ZeroTestProof>;

public:
    /*todo: откуда брать r?*/
    /*check: f = p -> f - p = 0*/
    static ptr_t forProver(CanonicPolynom &g, CanonicPolynom &p, TG_t tG, wGenerator_t wGenerator);
    static ptr_t forVerifier(commit_t comF, commit_t comQ, Y_t f_r, Y_t q_r);

    virtual bool check(wGenerator_t wGenerator = wGeneratorDefault) override;

private:
    ZeroTestProof() = default;

    TG_t m_tG;

    commit_t m_comF;
    PolynomSubstitutionProof m_fRproof;

    commit_t m_comQ;
    PolynomSubstitutionProof m_qRproof;

    X_t m_r;
    Y_t m_fR, m_qR;

    std::size_t m_witnessCount;
};

}

#endif // PROOF_H
