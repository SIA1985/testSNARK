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

    virtual bool check(G2 tG2, G2 g2) = 0;
};

class PlonkProof : public Proof
{
public:
    using ptr_t = std::shared_ptr<PlonkProof>;

    static ptr_t forProver(PartedCanonicPolynom &g, PartedCanonicPolynom &p, GPK_t &GPK, const witnesses_t &witness, Y_t u);
    static ptr_t forVerifier(commit_t comF, commit_t comQ, dot_t toProve, GPK_t GPK);

    virtual bool check(G2 tG2, G2 g2) override;

protected:
    virtual json_t toJson() const override;
    virtual bool fromJson(const json_t &json) override;

private:
    PlonkProof() = default;

    commit_t m_comF;
    commit_t m_comQ;

    dot_t m_toProve;

    GPK_t m_GPK;

    friend class ZeroTestProof;
    friend class ProverProof;
};

class PolynomSubstitutionProof : public Proof
{
public:
    using ptr_t = std::shared_ptr<PolynomSubstitutionProof>;

    static ptr_t forProver(CanonicPolynom &f, dot_t toProve, GPK_t GPK);
    static ptr_t forVerifier(commit_t comF, commit_t comQ, dot_t toProve, GPK_t GPK);

    virtual bool check(G2 tG2, G2 g2) override;

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

class ProverProof : public Proof
{
public:
    ProverProof() = default;

    ProverProof(const GlobalParams &gp, const values_t &input, value_t output);

    virtual bool check(G2 tG2, G2 g2) override;

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
    void currentOutput(PartedCanonicPolynom &t, value_t output, std::size_t lastWNum, GPK_t GPK);

};

}

#endif // PROOF_H
