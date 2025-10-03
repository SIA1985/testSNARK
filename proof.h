#ifndef PROOF_H
#define PROOF_H

#include <memory>

namespace snrk {

template <class Polynom>
class Proof
{
public:
    Proof() = default;
    virtual ~Proof() = default;

    virtual bool generate() = 0;

    virtual bool check() = 0;
};

template <class Polynom>
class ZeroTestProof : public Proof<Polynom>
{
public:
    ZeroTestProof(Polynom *f)
        : m_f{std::make_shared<Polynom>(f)}
    {
    }

private:
    std::shared_ptr<Polynom> m_f;
    std::shared_ptr<Polynom> m_g;
};

}

#endif // PROOF_H
