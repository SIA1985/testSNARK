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

ProverProof::ProverProof(const GlobalParams &gp, const values_t &input, value_t output)
{
    auto TParams = gp.PP().TParams;
    auto witnesses = gp.witnesses();

    m_GPK = gp.GPK();

    //todo: Зачем label в транскрипт?
    //Инициализация Транскрипта todo: (мб публичный ключ)
    Transcript tr("test");

    /*
   //1. Листья дерева (5 коммитментов в одном хеше) ---
    size_t numSegments = splA.size();
    std::vector<std::string> leafHashes;

    for(size_t i = 0; i < numSegments; ++i) {
        // Конкатенируем 7 точек в одну строку для хеширования
        std::string data = splA.getCommit(i).getStr()  + splB.getCommit(i).getStr()  +
                           splC.getCommit(i).getStr()  + splQ.getCommit(i).getStr()  +
                           splZ.getCommit(i).getStr()  + splW.getCommit(i).getStr()  +
                           splWT.getCommit(i).getStr();

        leafHashes.push_back(sha256(data));
    }

    MerkleTree tree(leafHashes);
    m_merkleRoot = tree.getRoot();
    tr.appendString("merkleRoot", m_merkleRoot);
    */

    //2. Получение точки раскрытия из Т
    auto r = tr("r");

    auto A = TParams.splittedT.left.toPartedCanonicPolynom();
    auto B = TParams.splittedT.right.toPartedCanonicPolynom();
    auto C = TParams.splittedT.result.toPartedCanonicPolynom();

    m_commitAr = A[r].commit(m_GPK);
    m_commitBr = B[r].commit(m_GPK);
    m_commitCr = C[r].commit(m_GPK);

    //3. Деление на Zero-полином (comQ -> T)
    auto p = correctGates(TParams.splittedT, gp.PP().SParams, gp.SWitnesses(), m_GPK);
    auto f = p - gp.PP().TParams.splittedT.result.toPartedCanonicPolynom();

    auto Z = ZeroWitnessPolynom(witnesses).toPartedCanonicPolynom();
    auto QG = f / Z;

    m_commitZr = Z[r].commit(m_GPK);
    m_commitQGr = QG[r].commit(m_GPK);

//    m_merklePath = tree.getProof(m_segmentIndex);
    //todo: фиксация в транскрипте корня Меркла?
    //todo: если r и r_next в разных сплайнах, то нужно их всех в Меркла


    //4. Полином-аккумулятор, для проверки перестановки
    auto beta = tr("beta");
    auto gamma = tr("gamma");

    auto WT = gp.PP().WParams.wt.toPartedCanonicPolynom();
    auto WI = gp.PP().WParams.wi.toPartedCanonicPolynom();

    dots_t WDots, WDotsShift1, NumDots, DenDots;
    Y_t currentW = 1;
    size_t n = gp.circutSize();

    WDots.push_back({0, currentW}); // W(0) = 1
    for (size_t i = 0; i < n; ++i) {
        value_t valA = A(i);
        value_t id   = WI(i);
        value_t sig  = WT(i);

        value_t num = valA + (beta * id) + gamma;
        value_t den = valA + (beta * sig) + gamma;

        currentW *= (num / den);

        WDots.push_back({(value_t)(i + 1), currentW});
        WDotsShift1.push_back({(value_t)(i), currentW});
        NumDots.push_back({});
        DenDots.push_back({});
    }
    WDotsShift1.push_back({(value_t)(n), 1});

    auto num = A + (WI * beta) + gamma;
    auto den = A + (WT * beta) + gamma;

    auto WShift1 = InterpolationPolynom(WDotsShift1).toPartedCanonicPolynom();
    auto W = InterpolationPolynom(WDots).toPartedCanonicPolynom();
    m_commitWr = W[r].commit(m_GPK);
    std::cout << WI.distance() << " " << WT.distance() << std::endl;

    auto Err = (WShift1 * num) - (W * den);
    auto QP = Err / Z;
    m_commitQPr = QP[r].commit(m_GPK);

    m_rA = A(r);
    m_rB = B(r);
    m_rC = C(r);
    m_rW = W(r);
    m_rWNext = W(r + 1);
    m_rZ = Z(r);
    m_rQG = QG(r);
    m_rQP = QP(r);

    m_rWT = WT(r);
    m_rWI = WI(r);

    tr.appendScalar("rA", m_rA);
    tr.appendScalar("rB", m_rB);
    tr.appendScalar("rC", m_rC);
    tr.appendScalar("rW", m_rW);
    tr.appendScalar("rWNext", m_rWNext);
    tr.appendScalar("rZ", m_rZ);
    tr.appendScalar("rQG", m_rQG);
    tr.appendScalar("rQP", m_rQP);

    //5. Получения доказательства pi
    auto v = tr("v");

    CanonicPolynom F;
    value_t currentV = 1;
    for(auto &poly : {(A[r] - m_rA), (B[r] - m_rB), (C[r] - m_rC),
                      (W[r] - m_rW), (Z[r] - m_rZ), (QG[r] - m_rQG),
                      (QP[r] - m_rQP)}) {
        F += poly * currentV;
        currentV *= v;
    }
    F += (W[r] - m_rWNext) * currentV;
    F -= F(r);

    auto z = CanonicPolynom({-r, 1});
    m_pi = (F / z).commit(m_GPK);
}

bool ProverProof::check(G2 tG2, G2 g2)
{
    //Инициализация Транскрипта todo: (мб публичный ключ)
    Transcript tr("test");

//    tr.appendString("merkleRoot", m_merkleRoot);
    auto r = tr("r");

    //Проверка пути Меркла
//    std::string leafData = m_CjA.getStr() + m_CjB.getStr() + m_CjC.getStr() + m_CjQ.getStr() + m_CjZ.getStr();
//    if (!verifyMerkle(m_merkleRoot, leafData, m_merklePath, m_segmentIndex)) return false;

    // Проверка уравнения связи W
    auto beta = tr("beta");
    auto gamma = tr("gamma");

    value_t left = m_rWNext * (m_rA + beta * m_rWI + gamma);
    value_t right = m_rW * (m_rA + beta * m_rWT + gamma);

    if (left - right != m_rQP * m_rZ) {
        std::cout << left - right << std::endl;
        std::cout << m_rZ << std::endl;
        std::cout << m_rQP << std::endl;
        return false;
    }

    tr.appendScalar("rA", m_rA);
    tr.appendScalar("rB", m_rB);
    tr.appendScalar("rC", m_rC);
    tr.appendScalar("rW", m_rW);
    tr.appendScalar("rWNext", m_rWNext);
    tr.appendScalar("rZ", m_rZ);
    tr.appendScalar("rQG", m_rQG);
    tr.appendScalar("rQP", m_rQP);

    auto v = tr("v");

    commit_t F;
    F.clear();
    value_t currentV = 1;
    for(auto &comm : {(m_commitAr - m_GPK.g1 * m_rA), (m_commitBr - m_GPK.g1 * m_rB),
                      (m_commitCr - m_GPK.g1 * m_rC), (m_commitWr - m_GPK.g1 * m_rW),
                      (m_commitZr - m_GPK.g1 * m_rZ), (m_commitQGr - m_GPK.g1 * m_rQG),
                      (m_commitQPr - m_GPK.g1 * m_rQP)}) {
        F += comm * currentV;
        currentV *= v;
    }
    F += (m_commitWr - m_GPK.g1 * m_rWNext) * currentV;

    value_t compensation = (m_rW - m_rWNext) * currentV;
    F -= (m_GPK.g1 * compensation);

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

PartedCanonicPolynom ProverProof::correctGates(const SplittedT_t &t, const GlobalParams::SParams_t SParams, const witnesses_t &ws, GPK_t GPK)
{
    const auto &left = t.left.toPartedCanonicPolynom();
    const auto &right = t.right.toPartedCanonicPolynom();

    auto funcF = PartedCanonicPolynom(PartedCanonicPolynom::map_t{});
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
