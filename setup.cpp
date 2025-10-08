#include "setup.h"

#include <algorithm>
#include <map>

namespace snrk {

W_t W_t::generate(const witnesses_t &witnesses, const Circut &circut)
{
    W_t w;
    std::map<value_t, std::shared_ptr<cond_t>> duplicates;
    auto insert = [&duplicates](value_t key, witness_t value)
    {
        if (!duplicates[key]) {
            duplicates[key] = std::make_shared<cond_t>();
        }

        duplicates[key]->insert(value);
    };

    auto cw = witnesses.cbegin();

    auto fillMap = [&cw, &duplicates, &insert](const values_t &row)
    {
        for(const auto &elment : row) {
            insert(elment, *cw++);
        }
    };

    fillMap(circut.m_inputX);
    fillMap(circut.m_inputW);

    for(const auto& gate : circut.m_gates) {
        insert(gate.m_input.a, *cw++);
        insert(gate.m_input.b, *cw++);
        insert(gate.m_output, *cw++);
    }

    for(const auto &[_, condition] : duplicates) {
        for(const auto &witness : *condition) {
            w.m_map[witness] = condition;
        }
    }

    return w;
}

W_t::cond_t W_t::operator()(witness_t w) const
{
    return *m_map.at(w);
}

GlobalParams::GlobalParams(const Circut &circut)
{
    witnesses_t witnesses;
    witnesses.resize(circut.degree());

    witness_t wGenerator = -circut.inputSize();
    std::generate(witnesses.begin(), witnesses.end(),
    [&wGenerator]() -> witness_t
    {
        return wGenerator++;
    });

    /*кадый в отдельный поток, т.к. читаем witness и circut*/
    generateT(witnesses, circut);
    generateS(circut);
    generateW(witnesses, circut);
}

GlobalParams::ProverParams GlobalParams::PP() const
{

}

GlobalParams::VerifierParams GlobalParams::VP() const
{

}

void GlobalParams::generateT(const witnesses_t &witnesses, const Circut &circut)
{
    typename T_t::dots_t dots;

    auto cw = witnesses.cbegin();
    auto fillMap = [&cw, &dots](const values_t &row)
    {
        for(const auto &elment : row) {
            dots.push_back({*cw++, elment});
        }
    };

    fillMap(circut.m_inputX);
    fillMap(circut.m_inputW);

    for(const auto& gate : circut.m_gates) {
        dots.push_back({*cw++, gate.m_input.a});
        dots.push_back({*cw++, gate.m_input.b});
        dots.push_back({*cw++, gate.m_output});
    }

    m_t = T_t::generate(dots);
}

void GlobalParams::generateS(const Circut &circut)
{
    typename S_t::dots_t dots;

    for(int i = 0; i < circut.size(); i++) {
        dots.push_back({i, circut.m_gates[i].m_type});
    }

    m_s = S_t::generate(dots);
}

void GlobalParams::generateW(const witnesses_t &witnesses, const Circut &circut)
{
    m_w = W_t::generate(witnesses, circut);
}

}
