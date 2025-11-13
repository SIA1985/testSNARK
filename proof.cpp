#include "proof.h"

#include <cmath>
#include <iostream>
#include <iomanip>

namespace snrk {

bool equal(const ValueType &a, const ValueType &b, double eps = 1e-9)
{
    ValueType c = a - b;
    if(c < 0) {
        c = -c;
    }

    std::cout << std::setprecision(20) << "diff " << a << " - " << b << " = " << c << std::endl;

    return c <= eps;
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

ZeroTestProof::ptr_t ZeroTestProof::forProver(PartedCanonicPolynom &g, PartedCanonicPolynom &p, TG_t tG, xs_t witness)
{
    auto ptr = ptr_t(new ZeroTestProof);

    ptr->m_tG = tG;

    auto f = g - p;

    ptr->m_comF = f.commit(tG);
    ptr->m_witness = witness;

    auto z = ZeroPolynom::generate(ptr->m_witness);

    auto q = f / z;

    ptr->m_comQ = q.commit(tG);

    /*todo: */
    ptr->m_r = 15;

    //todo: одинаковые значения!
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

    auto z = ZeroPolynom::generate(m_witness);

    ValueType b = m_qR * z(m_r);

    std::cout << std::setprecision(20) << m_fR << " " << m_qR << " " << z(m_r) << std::endl;
    return equal(m_fR, b);
}

}
