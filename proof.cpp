#include "proof.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include <nlohmann/json.hpp>
#include <thread>
#include <csignal>

namespace snrk {

bool equal(const value_t &a, const value_t &b, double eps = 1e-9)
{
    value_t c = a - b;
    if(c < value_t(0)) {
        c = -c;
    }

    return c <= eps;
}

X_t getR(const witnesses_t &witness, value_t t) {
    if (witness.size() == 0) {
        return 1; //todo
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    auto size = witness.size();

    witnesses_t::const_iterator randIt;
    std::uniform_int_distribution<decltype(size)> distrib(0, size - 1);
    int attemps = 3;

    do {
        randIt = witness.begin();
        std::advance(randIt, distrib(gen));

        if (--attemps == 0) {
            std::raise(SIGTERM);
        }
    }while(value_t(*randIt) == t);

    return *randIt;
}

PolynomSubstitutionProof::ptr_t PolynomSubstitutionProof::forProver(Polynom &f, dot_t toProve, TG_t tG)
{
    auto ptr = ptr_t(new PolynomSubstitutionProof);

    auto q = CustomPolynom([&f, u = toProve.x](X_t x) -> Y_t
    {
        return (f(x) - f(u)) / (x - u);
    });

    ptr->m_comF = f.commit(tG);
    ptr->m_comQ = q.commit(tG);
    ptr->m_toProve = toProve;
    ptr->m_tG = tG;

    return ptr;
}

PolynomSubstitutionProof::ptr_t PolynomSubstitutionProof::forVerifier(commit_t comF, commit_t comQ, dot_t toProve, TG_t tG)
{
    auto ptr = ptr_t(new PolynomSubstitutionProof);

    ptr->m_comF = comF;
    ptr->m_comQ = comQ;
    ptr->m_toProve = toProve;
    ptr->m_tG = tG;

    return ptr;
}

bool PolynomSubstitutionProof::check()
{
    value_t a = (m_tG.t - m_toProve.x) * m_comQ;
    value_t b = m_comF - m_toProve.y * m_tG.G;
    return equal(a, b);
}

json_t PolynomSubstitutionProof::toJson() const
{
    json_t json;

    ToJson(json, m_comF);
    ToJson(json, m_comQ);

    ToJson(json, m_toProve);
    ToJson(json, m_tG);

    return json;
}

bool PolynomSubstitutionProof::fromJson(const json_t &json)
{
    FromJson(json, m_comF);
    FromJson(json, m_comQ);
    FromJson(json, m_toProve);
    FromJson(json, m_tG);

    return true;
}

ZeroTestProof::ptr_t ZeroTestProof::forProver(PartedCanonicPolynom &g, PartedCanonicPolynom &p, TG_t tG, const witnesses_t &witness, witness_t wStep)
{
    assert(g.distance() == p.distance());


    auto ptr = ptr_t(new ZeroTestProof);

    ptr->m_tG = tG;

    auto f = g - p;
    ptr->m_comF = f.commit(tG);

    auto z = ZeroWitnessPolynom(witness).toPartedCanonicPolynom();

    auto q = f.mustDevide(z);

    ptr->m_comQ = q.commit(tG);

    /*todo: (hash % size(witness) + тут можно любое число!*/
    ptr->m_r = getR(witness, tG.t);

    auto rRange = z.atRange(ptr->m_r);

    witnesses_t rPratedWitnesses = genWitnesses(rRange.leftBound().get_ui(), PartedCanonicPolynom::Partition, wStep);
    ptr->m_rPartedWitnesses = rPratedWitnesses;

    ptr->m_fR = f(ptr->m_r);
    ptr->m_qR = q(ptr->m_r);

    ptr->m_fRproof = *PolynomSubstitutionProof::forProver(f, {ptr->m_r, ptr->m_fR}, ptr->m_tG);
    ptr->m_qRproof = *PolynomSubstitutionProof::forProver(q, {ptr->m_r, ptr->m_qR}, ptr->m_tG);

    return ptr;
}

bool ZeroTestProof::check()
{
    if (!m_fRproof.check()) {
        return false;
    }

    if (!m_qRproof.check()) {
        return false;
    }

    auto z = ZeroWitnessPolynom(m_rPartedWitnesses).toPartedCanonicPolynom();

    value_t b = m_qR * z(m_r);

    std::cout << std::setprecision(20) << m_fR << " " << m_qR << " " << z(m_r) << " : " << m_comQ << std::endl;
    return equal(m_fR, b);
}

json_t ZeroTestProof::toJson() const
{
    json_t json;

    ToJson(json, m_tG);

    ToJson(json, m_comF);
    ToJson(json, m_fRproof);

    ToJson(json, m_comQ);
    ToJson(json, m_qRproof);

    ToJson(json, m_r);
    ToJson(json, m_fR);
    ToJson(json, m_qR);

    ToJson(json, m_rPartedWitnesses);

    return json;
}

bool ZeroTestProof::fromJson(const json_t &json)
{
    FromJson(json, m_tG);

    FromJson(json, m_comF);
    FromJson(json, m_fRproof);

    FromJson(json, m_comQ);
    FromJson(json, m_qRproof);

    FromJson(json, m_r);
    FromJson(json, m_fR);
    FromJson(json, m_qR);

    FromJson(json, m_rPartedWitnesses);

    return true;
}


ProverProof::ProverProof(const GlobalParams &gp, const values_t &input, value_t output)
{
    auto TParams = gp.PP().TParams;
    auto witnesses = gp.witnesses();
    auto tG = gp.TG();
    auto tCanonic = TParams.t.toPartedCanonicPolynom();

    correctInputs(tCanonic, input, witnesses, tG);
    currentOutput(tCanonic, output, witnesses.size(), tG);

    #define r(a) std::ref(a)

    std::thread th1(&ProverProof::correctGates, this, r(TParams.splittedT), gp.PP().SParams, gp.SWitnesses(), tG);
    std::thread th2(&ProverProof::currentVars, this, gp.PP().WParams.wt, r(tCanonic), r(witnesses), tG);

    #undef r

    th1.join();
    th2.join();
}

bool ProverProof::check()
{
    if (!m_inputsProof.check()) {
        std::cout << "Некорректные входы!" << std::endl;
        return false;
    }

    if (!m_gatesProof.check()) {
        std::cout << "Некорректные переходы!" << std::endl;
        return false;
    }

    if (!m_varsProof.check()) {
        std::cout << "Некорректные переменные!" << std::endl;
        return false;
    }

    if (!m_outputProof.check()) {
        std::cout << "Некорректный выход!" << std::endl;
        return false;
    }

    return true;
}

json_t ProverProof::toJson() const
{
    json_t json;

    ToJson(json, m_inputsProof);
    ToJson(json, m_gatesProof);
    ToJson(json, m_varsProof);
    ToJson(json, m_outputProof);

    return json;
}

bool ProverProof::fromJson(const json_t &json)
{
    FromJson(json, m_inputsProof);
    FromJson(json, m_gatesProof);
    FromJson(json, m_varsProof);
    FromJson(json, m_outputProof);

    return true;
}

void ProverProof::correctInputs(const PartedCanonicPolynom &tCanonic, values_t inputs, const witnesses_t &ws, TG_t tG)
{
    dots_t inputsW;
    inputsW.reserve(inputs.size());
    for(std::size_t i = 0; i < inputs.size(); i++) {
        inputsW.push_back({ws[i], inputs[i]});
    }

    auto funcV = InterpolationPolynom(inputsW).toPartedCanonicPolynom();
    auto funcTCut = tCanonic.cut(funcV.distance());

    auto wStep = *(++ws.begin()) - ws.front();
    witnesses_t witness = genWitnesses(ws.front(), inputs.size(), wStep);
    m_inputsProof = *ZeroTestProof::forProver(funcTCut, funcV, tG, witness, wStep);
}

void ProverProof::correctGates(const SplittedT_t &t, const GlobalParams::SParams_t SParams, const witnesses_t &ws, TG_t tG)
{
    auto left = t.left.toPartedCanonicPolynom();
    auto right = t.right.toPartedCanonicPolynom();
    auto result = t.result.toPartedCanonicPolynom();
    //400ms - 10k свидетелей

    auto funcF = PartedCanonicPolynom(PartedCanonicPolynom::map{});
    for(const auto &[operation, dots] : SParams.opsFromS) {

        //
        auto isOperation = InterpolationPolynom(dots).toPartedCanonicPolynom();
        //150ms - 10k свидетелей!

        switch(operation) {
        case Sum: {
            funcF += (left + right) * isOperation;
            break;
        }
        case Product: {
            funcF += (left * right) * isOperation;
            break;
        }
        case Minus: {
            funcF += (left - right) * isOperation;
            break;
        }
            //todo:
        case Devide: {
//            funcF += left.mustDevide(right) * isOperation;
            break;
        }
        default:
            break;
        }
        //60ms - 10k свидетелей

    }
    //1100ms - 10k свидетелей

    auto wStep = *(++ws.begin()) - ws.front();
    //
    m_gatesProof = *ZeroTestProof::forProver(funcF, result, tG, ws, wStep);
    //150ms - 10k свидетелей
}

void ProverProof::currentVars(const WT_t &wt, PartedCanonicPolynom &tCanonic, const witnesses_t &ws, TG_t tG)
{
    auto wtCanonic = wt.toPartedCanonicPolynom();

    auto wStep = *(++ws.begin()) - ws.front();
    m_varsProof = *ZeroTestProof::forProver(tCanonic, wtCanonic, tG, ws, wStep);
}

void ProverProof::currentOutput(PartedCanonicPolynom &tCanonic, value_t output, std::size_t lastWNum, TG_t tG)
{
    auto outputDot = dot_t{lastWNum, output};

    m_outputProof = *PolynomSubstitutionProof::forProver(tCanonic, outputDot, tG);
}

}
