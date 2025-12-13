#ifndef PROOF_H
#define PROOF_H

#include "polynom.h"
#include "setup.h"

namespace snrk {

class Proof : public Jsonable
{
public:
    Proof() = default;
    virtual ~Proof() = default;

    virtual bool check() = 0;
};

class PolynomSubstitutionProof : public Proof
{
public:
    using commit_t = X_t;
    using ptr_t = std::shared_ptr<PolynomSubstitutionProof>;

    static ptr_t forProver(Polynom &f, dot_t toProve, TG_t tG);
    static ptr_t forVerifier(commit_t comF, commit_t comQ, dot_t toProve, TG_t tG);

    virtual bool check() override;

protected:
    virtual json_t toJson() const override;
    virtual bool fromJson(const json_t &json) override;

private:
    PolynomSubstitutionProof() = default;

    commit_t m_comF;
    commit_t m_comQ;

    dot_t m_toProve;

    TG_t m_tG;

    friend class ZeroTestProof;
    friend class ProverProof;
};

class ZeroTestProof : public Proof
{
public:
    using commit_t = X_t;
    using ptr_t = std::shared_ptr<ZeroTestProof>;

    static ptr_t forProver(PartedCanonicPolynom &g, PartedCanonicPolynom &p, TG_t tG, const witnesses_t &witness, witness_t wStep);
    static ptr_t forVerifier(commit_t comF, commit_t comQ, Y_t f_r, Y_t q_r);

    virtual bool check() override;

protected:
    virtual json_t toJson() const override;
    virtual bool fromJson(const json_t &json) override;

private:
    ZeroTestProof() = default;

    TG_t m_tG;

    commit_t m_comF;
    PolynomSubstitutionProof m_fRproof;

    commit_t m_comQ;
    PolynomSubstitutionProof m_qRproof;

    X_t m_r;
    Y_t m_fR, m_qR;

    witnesses_t m_rPartedWitnesses;

    friend class ProverProof;
};

class ProverProof : public Proof
{
public:
    ProverProof(const GlobalParams &gp);

    virtual bool check() override;

protected:
    virtual json_t toJson() const override;
    virtual bool fromJson(const json_t &json) override;

private:
    ZeroTestProof m_inputsProof;
    ZeroTestProof m_gatesProof;
    ZeroTestProof m_varsProof;
    PolynomSubstitutionProof m_outputProof;
};

}

#endif // PROOF_H
