#include "proof.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include <nlohmann/json.hpp>
#include <thread>
#include <csignal>

namespace snrk {

Transcript::Transcript(const std::string &init)
    : m_state{init}
{
}

void Transcript::appendPoint(const std::string &label, const mcl::G1 &point)
{
    m_state += label;
    char buf[1024];
    size_t n = point.serialize(buf, sizeof(buf));
    updateState({buf, n});
}

void Transcript::appendScalar(const std::string &label, const mcl::Fr &scalar)
{
    m_state += label;
    char buf[1024];
    size_t n = scalar.serialize(buf, sizeof(buf));
    updateState({buf, n});
}

value_t Transcript::operator()(const std::string &label)
{
    m_state += label;

    value_t challenge;
    challenge.setHashOf(m_state);

    updateState(m_state);

    return challenge;
}

void Transcript::updateState(const std::string &data)
{
    m_state.append(data);

    mcl::bn::Fr hash_value;
    hash_value.setHashOf(m_state);
    m_state = hash_value.getStr(16);
}


PlonkProof::ptr_t PlonkProof::forProver(PartedCanonicPolynom &g, PartedCanonicPolynom &p, GPK_t &GPK, const witnesses_t &witness, X_t u)
{
    assert(g.distance() == p.distance());


    auto ptr = ptr_t(new PlonkProof);

    ptr->m_GPK = GPK;

    auto f = g - p;
    ptr->m_comF = f.commit(GPK);

    auto z = ZeroWitnessPolynom(witness).toPartedCanonicPolynom();

    auto q = f / z;
    ptr->m_comQ = q.commit(GPK);

    ptr->m_toProve.x = u;
    ptr->m_toProve.y = q(u);

    return ptr;
}

PlonkProof::ptr_t PlonkProof::forVerifier(commit_t comF, commit_t comQ, dot_t toProve, GPK_t GPK)
{
    auto ptr = ptr_t(new PlonkProof);

    ptr->m_comF = comF;
    ptr->m_comQ = comQ;
    ptr->m_toProve = toProve;
    ptr->m_GPK = GPK;

    return ptr;
}

bool PlonkProof::check(G2 tG2, G2 g2)
{
    G2 zG2, forE1;
    G2::mul(zG2, g2, m_toProve.x);
    G2::sub(forE1, tG2, zG2);

    GT e1;
    mcl::pairing(e1, m_comQ, forE1);

    G1 yG, forE2;
    G1::mul(yG, m_GPK.g, m_toProve.y);
    G1::sub(forE2, m_comF, yG);

    GT e2;
    mcl::pairing(e2, forE2, g2);

    return e1 == e2;
}

json_t PlonkProof::toJson() const
{
    //todo:
    return {};
}

bool PlonkProof::fromJson(const json_t &json)
{
    //todo:
    return false;
}

ProverProof::ProverProof(const GlobalParams &gp, const values_t &input, value_t output)
{
    auto TParams = gp.PP().TParams;
    auto witnesses = gp.witnesses();

    m_GPK = gp.GPK();

    //todo: Зачем label в транскрипт?

    //1. Инициализация Транскрипта todo: (мб публичный ключ)
    Transcript tr("test");

    //2. Получаем обязательства splittedT и добавляем в T

    auto left = TParams.splittedT.left.toPartedCanonicPolynom();
    m_commitA = left.commit(m_GPK);

    auto right = TParams.splittedT.right.toPartedCanonicPolynom();
    m_commitB = right.commit(m_GPK);

    auto result = TParams.splittedT.result.toPartedCanonicPolynom();
    m_commitC = result.commit(m_GPK);

    tr.appendPoint("commitA", m_commitA);
    tr.appendPoint("commitB", m_commitB);
    tr.appendPoint("commitC", m_commitC);


    //3. Какой-то полином-аккумулятор, что у меня является W(x) и WT(x)

    //4. Деление на Zero-полином (comQ -> T)
    auto p = correctGates(TParams.splittedT, gp.PP().SParams, gp.SWitnesses(), m_GPK);
    auto f = p - gp.PP().TParams.splittedT.result.toPartedCanonicPolynom();
//    ptr->m_comF = f.commit(GPK);

    auto z = ZeroWitnessPolynom(witnesses).toPartedCanonicPolynom();
    auto q = f / z;

    m_commitQ = q.commit(m_GPK);
    tr.appendPoint("commitQ", m_commitQ);

    //5. Получение точки раскрытия из Т
    m_r = tr("r");

    //6. Получение значений всех полиномов выше (a, b, c, w, q) (для последующей проверки в SubstitutionProof,
    //   точнее он будет фигурировать в проверке в "разобранном" виде)
    m_rA = TParams.splittedT.left(m_r);
    m_rB = TParams.splittedT.right(m_r);
    m_rC = TParams.splittedT.result(m_r);
//    auto Rw;
    m_rQ = q(m_r);

    m_commitAQ = getCommitQ(left[m_r], m_r, m_GPK);
    m_commitBQ = getCommitQ(right[m_r], m_r, m_GPK);
    m_commitCQ = getCommitQ(result[m_r], m_r, m_GPK);
}

bool ProverProof::check(G2 tG2, G2 g2)
{
    auto result = true;

    auto aProof = PlonkProof::forVerifier(m_commitA, m_commitAQ, {m_r, m_rA}, m_GPK);
    auto bProof = PlonkProof::forVerifier(m_commitB, m_commitBQ, {m_r, m_rB}, m_GPK);
    auto cProof = PlonkProof::forVerifier(m_commitC, m_commitCQ, {m_r, m_rC}, m_GPK);
    //tg2, g2 для проверки pi-пруфа

    for(auto &proof : {aProof, bProof, cProof}) {
        result = result && proof->check(tG2, g2);
    }

    return result;
}

json_t ProverProof::toJson() const
{
    //todo:
    json_t json;

    return json;
}

bool ProverProof::fromJson(const json_t &json)
{
    //todo:
    return false;
}

commit_t ProverProof::getCommitQ(CanonicPolynom f, X_t toProve, GPK_t GPK) const
{
    auto u = toProve;
    auto v = f(u);
    auto xSubU = CanonicPolynom({-u, 1});

    auto q = (f - CanonicPolynom({v, 0})) / xSubU;

    return q.commit(GPK);
}

PartedCanonicPolynom ProverProof::correctGates(const SplittedT_t &t, const GlobalParams::SParams_t SParams, const witnesses_t &ws, GPK_t GPK)
{
    const auto &left = t.left.toPartedCanonicPolynom();
    const auto &right = t.right.toPartedCanonicPolynom();

    auto funcF = PartedCanonicPolynom(PartedCanonicPolynom::map{});
    for(const auto &[operation, dots] : SParams.opsFromS) {

        auto isOperation = InterpolationPolynom(dots).toPartedCanonicPolynom();

        switch(operation) {
        case Sum: {
            funcF += (left + right) * isOperation;
            break;
        }
        case Product: {
            funcF += (left * right) * isOperation;
            break;
        }
        default:
            assert(false);
        }

    }

    return funcF;
}

}
