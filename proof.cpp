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

void Transcript::appendHash(const std::string &label, const hash_t &hashData)
{
    m_state += label;

    updateState(hashData);
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
    m_state = hash(m_state + data);
}


PlonkProof::ptr_t PlonkProof::forProver(SplinePolynom &g, SplinePolynom &p, GPK_t &GPK, const witnesses_t &witness, X_t u)
{
    assert(g.distance() == p.distance());


    auto ptr = ptr_t(new PlonkProof);

    ptr->m_GPK = GPK;

    auto f = g - p;
    //todo:
//    ptr->m_comF = f.commit(GPK);

    auto z = ZeroWitnessPolynom(witness).toSplinePolynom();

    auto q = f / z;
    //todo:
//    ptr->m_comQ = q.commit(GPK);

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
    G1::mul(yG, m_GPK.g1, m_toProve.y);
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

ProverProof::ProverProof(const GlobalParams &gp)
{
    auto TParams = gp.PP().TParams;
    auto witnesses = gp.witnesses();

    m_GPK = gp.GPK();

    //todo: Зачем label в транскрипт?
    //1. Инициализация Транскрипта todo: (мб публичный ключ)
    Transcript tr("test");

    //2. Построение всех сплайнов
    auto Z = ZeroWitnessPolynom(witnesses).toSplinePolynom();

    auto T = TParams.t.toSplinePolynom();

    auto p = correctGates(TParams.splittedT, gp.PP().SParams);
    auto f = p - gp.PP().TParams.splittedT.result.toSplinePolynom();

    auto QG = f / Z;

    auto beta = tr("beta");
    auto gamma = tr("gamma");

    auto WT = gp.PP().WParams.wt.toSplinePolynom();
    auto WI = gp.PP().WParams.wi.toSplinePolynom();

    auto num = T + (WT * beta) + gamma;
    auto den = T + (WI * beta) + gamma;

    auto [W, WShift1] = correctPermulations(witnesses, num, den);

    auto Err = (WShift1 * den) - (W * num);
    auto QP = Err / Z;


    //3. Дерево Меркла для доказательства монолитности сплайнов
     size_t numSegments = Z.segmentsCounts();
     hashes_t leafHashes;

     for(size_t i = 0; i < numSegments; ++i) {
         std::string data;
         for(const auto &poly : {T, W, WShift1, Z, QG, QP, WT, WI}) {
             data += poly.commit(m_GPK, i).getStr();
         }

         leafHashes.push_back(hash(data));
     }

     MerkleTree tree(leafHashes);
     m_merkleRoot = tree.root();
     tr.appendHash("merkleRoot", m_merkleRoot);

    //4. Получение точки раскрытия из Т
    auto r = tr("r");

    std::string toHash;
    for(const auto &poly : {T[r], W[r], WShift1[r], Z[r],
                            QG[r], QP[r], WT[r], WI[r]}) {
        toHash += poly.commit(m_GPK).getStr();
    }
    m_merkleLeaf = hash(toHash);

    m_merklePath = tree.path(m_merkleLeaf).value();
    //todo: если r и r_next в разных сплайнах, то нужно их всех в Меркла?

    m_commitTr = T[r].commit(m_GPK);
    m_commitZr = Z[r].commit(m_GPK);
    m_commitQGr = QG[r].commit(m_GPK);
    m_commitWr = W[r].commit(m_GPK);
    m_commitWNextr = WShift1[r].commit(m_GPK);
    m_commitQPr = QP[r].commit(m_GPK);

    m_rT = T(r);
    m_rW = W(r);
    m_rWNext = WShift1(r);
    m_rZ = Z(r);
    m_rQG = QG(r);
    m_rQP = QP(r);

    m_rWT = WT(r);
    m_rWI = WI(r);

    tr.appendScalar("rT", m_rT);
    tr.appendScalar("rW", m_rW);
    tr.appendScalar("rWNext", m_rWNext);
    tr.appendScalar("rZ", m_rZ);
    tr.appendScalar("rQG", m_rQG);
    tr.appendScalar("rQP", m_rQP);


    //5. Получения доказательства pi
    auto v = tr("v");

    CanonicPolynom F;
    value_t currentV = 1;
    for(auto &poly : {(T[r] - m_rT), (W[r] - m_rW), (WShift1[r] - m_rWNext),
                      (Z[r] - m_rZ), (QG[r] - m_rQG), (QP[r] - m_rQP)}) {
        F += poly * currentV;
        currentV *= v;
    }

    auto z = CanonicPolynom({-r, 1});
    m_pi = (F / z).commit(m_GPK);
}

bool ProverProof::check(G2 tG2, G2 g2)
{
    //Инициализация Транскрипта todo: (мб публичный ключ)
    Transcript tr("test");

    auto beta = tr("beta");
    auto gamma = tr("gamma");

    value_t numR = m_rT + (beta * m_rWT) + gamma;
    value_t denR = m_rT + (beta * m_rWI) + gamma;

    value_t errorW = (m_rWNext * denR) - (m_rW * numR);

    if (errorW != m_rQP * m_rZ) {
        return false;
    }

    tr.appendHash("merkleRoot", m_merkleRoot);

    auto r = tr("r");

    //Проверка пути Меркла
//    std::string leafData = m_CjA.getStr() + m_CjB.getStr() + m_CjC.getStr() + m_CjQ.getStr() + m_CjZ.getStr();
//    if (!verifyMerkle(m_merkleRoot, leafData, m_merklePath, m_segmentIndex)) return false;


    tr.appendScalar("rT", m_rT);
    tr.appendScalar("rW", m_rW);
    tr.appendScalar("rWNext", m_rWNext);
    tr.appendScalar("rZ", m_rZ);
    tr.appendScalar("rQG", m_rQG);
    tr.appendScalar("rQP", m_rQP);

    auto v = tr("v");

    commit_t F;
    F.clear();
    value_t currentV = 1;
    for(auto &comm : {(m_commitTr - m_GPK.g1 * m_rT), (m_commitWr - m_GPK.g1 * m_rW), (m_commitWNextr - m_GPK.g1 * m_rWNext),
                      (m_commitZr - m_GPK.g1 * m_rZ), (m_commitQGr - m_GPK.g1 * m_rQG),(m_commitQPr - m_GPK.g1 * m_rQP)}) {
        F += comm * currentV;
        currentV *= v;
    }

    GT e1, e2;
    mcl::pairing(e1, m_pi, tG2 - g2 * r);
    mcl::pairing(e2, F, g2);

    return e1 == e2;
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

SplinePolynom ProverProof::correctGates(const SplittedT_t &t, const GlobalParams::SParams_t SParams)
{
    const auto &left = t.left.toSplinePolynom();
    const auto &right = t.right.toSplinePolynom();

    auto funcF = SplinePolynom(SplinePolynom::map_t{});
    for(const auto &[operation, dots] : SParams.opsFromS) {

        auto isOperation = InterpolationPolynom(dots).toSplinePolynom();

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

ProverProof::WResult_t ProverProof::correctPermulations(const witnesses_t &witnesses, SplinePolynom &num, SplinePolynom &den)
{
    dots_t WDots, WDotsShift1;
    Y_t currN, currD, currW = 1;

    WDots.push_back({witnesses.front(), 1});
    for (auto it = witnesses.begin(); it != std::prev(witnesses.end()); it++) {
        currN = num(*it);
        currD = den(*it);

        currW *= currN / currD;
        WDots.push_back({*(it + 1), currW});
        WDotsShift1.push_back({*it, currW});
    }

    WDotsShift1.push_back({witnesses.back(), 1});

    return {InterpolationPolynom(WDots).toSplinePolynom(),
            InterpolationPolynom(WDotsShift1).toSplinePolynom()};
}

}
