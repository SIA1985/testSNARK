#ifndef PROOF_H
#define PROOF_H

#include "polynom.h"
#include "setup.h"

namespace snrk {

//todo: использовать другой хэш
class Transcript {
public:
   Transcript(const std::string &init);

   void appendPoint(const std::string &label, const mcl::bn::G1 &point);

   void appendScalar(const std::string &label, const mcl::bn::Fr &scalar);

   value_t operator()(const std::string &label);

private:
   void updateState(const std::string &data);

   std::string m_state;
};

class Proof : public Jsonable
{
public:
    Proof() = default;
    virtual ~Proof() = default;

    virtual bool check(G2 tG2, G2 g2) = 0;
};

/*Нуль-тест + случайная точка*/
class PlonkProof : public Proof
{
public:
    using ptr_t = std::shared_ptr<PlonkProof>;

    static ptr_t forProver(SplinePolynom &g, SplinePolynom &p, GPK_t &GPK, const witnesses_t &witness, X_t u);
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

    friend class ProverProof;
};

class ProverProof : public Proof
{
    using WResult_t = struct{SplinePolynom W; SplinePolynom WShift1;};
public:
    ProverProof() = default;

    ProverProof(const GlobalParams &gp);

    virtual bool check(G2 tG2, G2 g2) override;

protected:
    virtual json_t toJson() const override;
    virtual bool fromJson(const json_t &json) override;

private:
    //Подготовка увеличивается согласно О(n^2)
    SplinePolynom correctGates(const SplittedT_t &t, const GlobalParams::SParams_t SParams);
    WResult_t correctPermulations(const witnesses_t &witnesses, SplinePolynom &num, SplinePolynom &den);

    commit_t m_commitTr, m_commitWr, m_commitWNextr, m_commitZr, m_commitQGr, m_commitQPr;
    commit_t m_pi;

    Y_t m_rT, m_rW, m_rWNext, m_rZ, m_rQG, m_rQP;
    Y_t m_rWT, m_rWI;

    hash_t m_MerkleRoot;
    hashes_t m_MerklePath;

    //todo: по идее не нужен здесь, ибо глобальные параметры
    //мб добавить g2 тоже в GPK?
    GPK_t m_GPK;
};

}

#endif // PROOF_H
