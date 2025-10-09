#include "proof.h"

snrk::PolynomSubstitutionProof::PolynomSubstitutionProof(commit_t comF, commit_t comQ, dot_t toProve)
    : m_comF{comF}
    , m_comQ{comQ}
    , m_toProve{toProve}
{

}

void snrk::PolynomSubstitutionProof::setGp(X_t t, int G)
{
    m_t = t;
    m_G = G;
}

bool snrk::PolynomSubstitutionProof::check()
{
    /*todo: сравнение double*/
    return (m_t - m_toProve.x) * m_comQ == m_comF - m_toProve.y * m_G;
}
