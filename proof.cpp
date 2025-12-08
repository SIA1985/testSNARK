#include "proof.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

namespace snrk {

bool equal(const ValueType &a, const ValueType &b, double eps = 1e-9)
{
    ValueType c = a - b;
    if(c < 0) {
        c = -c;
    }

    return c <= eps;
}

X_t getR(const xs_t &witness) {
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
    ValueType a = (m_tG.t - m_toProve.x) * m_comQ;
    ValueType b = m_comF - m_toProve.y * m_tG.G;
    return equal(a, b);
}

#define printDur(text, end, start)     std::cout << text << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

ZeroTestProof::ptr_t ZeroTestProof::forProver(PartedCanonicPolynom &g, PartedCanonicPolynom &p, TG_t tG, xs_t witness)
{
    assert(g.distance() == p.distance());


    auto ptr = ptr_t(new ZeroTestProof);

    ptr->m_tG = tG;

    ptr->m_witness = witness;

    auto f = g - p;
    ptr->m_comF = f.commit(tG);

    //900ms при 30к свидетелей, 100ms - 10k
    auto z = ZeroPolynom(ptr->m_witness).toPartedCanonicPolynom();

//    for(auto w : ptr->m_witness) {
//        std::cout << w << " : " << g(w) << " - " << p(w) << std::endl;
//    }


    auto q = f.mustDevide(z);

    ptr->m_comQ = q.commit(tG);

    /*todo: (hash % size(witness) + тут можно любое число!*/
    ptr->m_r = getR(ptr->m_witness);

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

    //todo: оптимизация - создание zero только для промежутка, где есть m_r
    auto z = ZeroPolynom(m_witness).toPartedCanonicPolynom();

    ValueType b = m_qR * z(m_r);

//    std::cout << std::setprecision(20) << m_fR << " " << m_qR << " " << z(m_r) << std::endl;
    return equal(m_fR, b);
}

}
