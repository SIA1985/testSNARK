#include "proof.h"

#include <cmath>

namespace snrk {

bool equal(double a, double b, double eps = 1e-9)
{
    return std::fabs(a - b) <= eps;
}

xs_t genWitnessXs(std::size_t count)
{
    xs_t witnesses;

    for(std::size_t i = 1; i <= count; i++) {
        witnesses.insert(i);
    }

    return witnesses;
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
    auto a = (m_tG.t - m_toProve.x) * m_comQ;
    auto b = m_comF - m_toProve.y * m_tG.G;
    return equal(a, b);
}

ZeroTestProof::ptr_t ZeroTestProof::forProver(CanonicPolynom &g, CanonicPolynom &p, TG_t tG)
{
    auto ptr = ptr_t(new ZeroTestProof);

    ptr->m_tG = tG;

    auto f = g - p;

    ptr->m_comF = f.commit(tG);
    ptr->m_witnessCount = p.degree() + 1;

    auto z = ZeroPolynom::generate(genWitnessXs(ptr->m_witnessCount));

    auto q = f / z;

    ptr->m_comQ = q.commit(tG);

    /*todo: */
    ptr->m_r = 15;
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

    auto z = ZeroPolynom::generate(genWitnessXs(m_witnessCount));

    auto a = m_fR / z(m_r);
    auto b = m_qR;
    return equal(a, b);
}

}
