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

    virtual bool check(mcl::G2 tG2, mcl::G2 G2) = 0;
};

class PolynomSubstitutionProof : public Proof
{
public:
    using ptr_t = std::shared_ptr<PolynomSubstitutionProof>;

    static ptr_t forProver(CanonicPolynom &f, dot_t toProve, GPK_t GPK);
    static ptr_t forVerifier(commit_t comF, commit_t comQ, dot_t toProve, GPK_t GPK);

    virtual bool check(mcl::G2 tG2, mcl::G2 G2) override;

protected:
    virtual json_t toJson() const override;
    virtual bool fromJson(const json_t &json) override;

private:
    PolynomSubstitutionProof() = default;

    commit_t m_comF;
    commit_t m_comQ;

    dot_t m_toProve;

    GPK_t m_GPK;

    friend class ZeroTestProof;
    friend class ProverProof;
};

class ZeroTestProof : public Proof
{
public:
    using ptr_t = std::shared_ptr<ZeroTestProof>;

    static ptr_t forProver(PartedCanonicPolynom &g, PartedCanonicPolynom &p, GPK_t GPK, const witnesses_t &witness, witness_t wStep);
    static ptr_t forVerifier(commit_t comF, commit_t comQ, Y_t f_r, Y_t q_r);

    virtual bool check(mcl::G2 tG2, mcl::G2 G2) override;

protected:
    virtual json_t toJson() const override;
    virtual bool fromJson(const json_t &json) override;

private:
    ZeroTestProof() = default;

    GPK_t m_GPK;

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
    ProverProof() = default;

    ProverProof(const GlobalParams &gp, const values_t &input, value_t output);

    virtual bool check(mcl::G2 tG2, mcl::G2 G2) override;

protected:
    virtual json_t toJson() const override;
    virtual bool fromJson(const json_t &json) override;

private:
    //Подготовка увеличивается согласно О(n^2)
    void correctInputs(const PartedCanonicPolynom &tCanonic, values_t inputs, const witnesses_t &ws, GPK_t GPK);
    //Подготовка увеличивается согласно О(n^2)
    void correctGates(const SplittedT_t &t, const GlobalParams::SParams_t SParams, const witnesses_t &ws, GPK_t GPK);
    //мб дело в том, что надо проверять не t, а такое t, что выводит адреса
    void currentVars(const WT_t &wt,  PartedCanonicPolynom &tCanonic, const witnesses_t &ws, GPK_t GPK);
//    void currentOutput(PartedCanonicPolynom &t, value_t output, std::size_t lastWNum, GPK_t GPK);



    ZeroTestProof m_inputsProof;
    ZeroTestProof m_gatesProof;
    ZeroTestProof m_varsProof;
    PolynomSubstitutionProof m_outputProof;
};

}

#endif // PROOF_H
