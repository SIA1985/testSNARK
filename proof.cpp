#include "proof.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>

namespace snrk {

bool equal(const ValueType &a, const ValueType &b, double eps = 1e-9)
{
    ValueType c = a - b;
    if(c < 0) {
        c = -c;
    }

//    std::cout << std::setprecision(20) << "diff " << a << " - " << b << " = " << c << std::endl;

    return c <= eps;
}

X_t getR(const xs_t &witness) {
    if (witness.size() == 0) {
        return 0; //todo
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

    auto q = snrk::CustomPolynom::generate([&f, u = toProve.x](snrk::X_t x) -> snrk::Y_t
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

//todo: если заработает witness: xs_t -> witnesses_t
ZeroTestProof::ptr_t ZeroTestProof::forProver(PartedCanonicPolynom &g, PartedCanonicPolynom &p, TG_t tG, xs_t witness)
{
    auto ptr = ptr_t(new ZeroTestProof);

    ptr->m_tG = tG;

//    witness.erase(std::prev(witness.end()));

//    auto fDots = (g - p).dots(witness);
//    auto newWitnesses = genWitnesses(wStart, witness.size(), true);

//    witness.clear();
//    for(std::size_t i = 0; i < newWitnesses.size(); i++) {
//        witness.insert(newWitnesses[i]);
//        fDots[i].x = newWitnesses[i];
//    }


//    auto f = InterpolationPolynom::generate(fDots).toCanonicPolynom();

    auto f = g - p;

    ptr->m_comF = f.commit(tG);
    ptr->m_witness = witness;

    auto z = ZeroPolynom::generate(ptr->m_witness).toPartedCanonicPolynom();
//    auto z = CanonicPolynom::generate(CanonicPolynom::coefsFromRoots(ptr->m_witness));
//    for(auto w : ptr->m_witness) {
//        std::cout << w << " " <<  f(w + 0.1) << " / " << z(w) << std::endl;
//    }

    auto q = f / z;

    ptr->m_comQ = q.commit(tG);

    /*todo: (hash % size(witness)*/
    ptr->m_r = getR(ptr->m_witness);
    std::cout << "R: " << ptr->m_r << " " << f( ptr->m_r) << " / " << z( ptr->m_r) << " =?= " << q( ptr->m_r)<< std::endl;

//    for(auto w : ptr->m_witness) {
//        std::cout << w << " : " << f(w) << " / " << z(w) << " =?= " << q(w) << std::endl;
//    }

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

    auto z = ZeroPolynom::generate(m_witness).toPartedCanonicPolynom();
//    auto z = CanonicPolynom::generate(CanonicPolynom::coefsFromRoots(m_witness));

    ValueType b = m_qR * z(m_r);

    std::cout << std::setprecision(20) << m_fR << " " << m_qR << " " << z(m_r) << std::endl;
    return equal(m_fR, b);
}

}
