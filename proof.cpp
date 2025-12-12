#include "proof.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include "nlohmann/json.hpp"

namespace snrk {

bool equal(const value_t &a, const value_t &b, double eps = 1e-9)
{
    value_t c = a - b;
    if(c < value_t(0)) {
        c = -c;
    }

    return c <= eps;
}

X_t getR(const witnesses_t &witness) {
    if (witness.size() == 0) {
        return 1; //todo
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    auto size = witness.size();

    std::uniform_int_distribution<decltype(size)> distrib(0, size - 1);
    auto random_offset = distrib(gen);

    auto randIt = witness.begin();
    std::advance(randIt, random_offset);

    return *randIt;
}

mp_exp_t e = -10;
#define name(a) #a

#define containsJson(json, a) json.contains(name(a))

#define varToJson(json, var)        json[name(var)] = var
#define varFromJson(json, var)      var = json[name(var)].get<decltype(var)>()

#define valueToJson(json, value)    json[name(value)] = value.get_str(e)
#define valueFromJson(json, value)  value = (std::string)json[name(value)]

#define dotToJson(json, dot) {json_t _j; _j["x"] = dot.x.get_str(e); _j["y"] = dot.y.get_str(e); json[name(dot)] = _j;}
#define dotFromJson(json, dot) {dot.x = (std::string)json[name(dot)]["x"]; dot.y = (std::string)json[name(dot)]["y"];}

#define tgToJson(json, tg) {json_t _j; _j["t"] = tg.t.get_str(e); _j["G"] = tg.G; json[name(tg)] = _j;}
#define tgFromJson(json, tg) {tg.t = (std::string)json[name(tg)]["t"]; tg.G = json[name(tg)]["G"];}

#define jsonableToJson(json, j) json[name(j)] = j.toJson()
#define jsonableFromJson(json, j) j.fromJson(json[name(j)])

Proof::Proof(json_t json)
{
    fromJson(json);
}

PolynomSubstitutionProof::ptr_t PolynomSubstitutionProof::forProver(Polynom &f, dot_t toProve, TG_t tG)
{
    auto ptr = ptr_t(new PolynomSubstitutionProof);

    auto q = snrk::CustomPolynom([&f, u = toProve.x](snrk::X_t x) -> snrk::Y_t
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

    valueToJson(json, m_comF);
    valueToJson(json, m_comQ);

    dotToJson(json, m_toProve);
    tgToJson(json, m_tG);

    return json;
}

bool PolynomSubstitutionProof::fromJson(const json_t &json)
{
    if (!containsJson(json, m_comF)) {
        return false;
    }
    valueFromJson(json, m_comF);

    if (!containsJson(json, m_comQ)) {
        return false;
    }
    valueFromJson(json, m_comQ);

    if (!containsJson(json, m_toProve)) {
        return false;
    }
    dotFromJson(json, m_toProve);

    if (!containsJson(json, m_tG)) {
        return false;
    }
    tgFromJson(json, m_tG);

    return true;
}

ZeroTestProof::ptr_t ZeroTestProof::forProver(PartedCanonicPolynom &g, PartedCanonicPolynom &p, TG_t tG, const witnesses_t &witness, witness_t wStep)
{
    assert(g.distance() == p.distance());


    auto ptr = ptr_t(new ZeroTestProof);

    ptr->m_tG = tG;

    auto f = g - p;
    ptr->m_comF = f.commit(tG);

    //900ms при 30к свидетелей, 100ms - 10k
    auto z = ZeroWitnessPolynom(witness).toPartedCanonicPolynom();

//    for(auto w : witness) {
//        std::cout << w << " : " << g(w) << " - " << p(w) << " = " << f(w) << std::endl;
//    }

    auto q = f.mustDevide(z);

    ptr->m_comQ = q.commit(tG);

    /*todo: (hash % size(witness) + тут можно любое число!*/
    ptr->m_r = getR(witness);

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

    tgToJson(json, m_tG);

    valueToJson(json, m_comF);
    jsonableToJson(json, m_fRproof);

    valueToJson(json, m_comQ);
    jsonableToJson(json, m_qRproof);

    valueToJson(json, m_r);
    valueToJson(json, m_fR);
    valueToJson(json, m_qR);

    varToJson(json, m_rPartedWitnesses);

    return json;
}

bool ZeroTestProof::fromJson(const json_t &json)
{
    if (!containsJson(json, m_tG)) {
        return false;
    }
    tgFromJson(json, m_tG);

    if (!containsJson(json, m_comF)) {
        return false;
    }
    valueFromJson(json, m_comF);

    if (!containsJson(json, m_fRproof)) {
        return false;
    }
    jsonableFromJson(json, m_fRproof);

    if (!containsJson(json, m_comQ)) {
        return false;
    }
    valueFromJson(json, m_comQ);

    if (!containsJson(json, m_qRproof)) {
        return false;
    }
    jsonableFromJson(json, m_qRproof);

    if (!containsJson(json, m_r)) {
        return false;
    }
    valueFromJson(json, m_r);

    if (!containsJson(json, m_fR)) {
        return false;
    }
    valueFromJson(json, m_fR);

    if (!containsJson(json, m_qR)) {
        return false;
    }
    valueFromJson(json, m_qR);

    if (!containsJson(json, m_rPartedWitnesses)) {
        return false;
    }
    varFromJson(json, m_rPartedWitnesses);

    return true;
}


ProverProof::ProverProof(const GlobalParams &gp)
{

}

bool ProverProof::check()
{
    return false;
}

json_t ProverProof::toJson() const
{
    json_t json;

    jsonableToJson(json, m_inputsProof);
    jsonableToJson(json, m_gatesProof);
    jsonableToJson(json, m_varsProof);
    jsonableToJson(json, m_outputProof);

    return json;
}

bool ProverProof::fromJson(const json_t &json)
{
    if (!containsJson(json, m_inputsProof)) {
        return false;
    }
    jsonableFromJson(json, m_inputsProof);

    if (!containsJson(json, m_gatesProof)) {
        return false;
    }
    jsonableFromJson(json, m_gatesProof);

    if (!containsJson(json, m_varsProof)) {
        return false;
    }
    jsonableFromJson(json, m_varsProof);

    if (!containsJson(json, m_outputProof)) {
        return false;
    }
    jsonableFromJson(json, m_outputProof);

    return true;
}

}
