#include "setup.h"
#include "operation.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <thread>

namespace snrk {

witnesses_t genWitnesses(witness_t start, std::size_t count, bool Chebishev)
{
    witnesses_t witnesses;
    for(std::size_t i = 1; i <= count; i++) {
        witnesses.push_back(
            Chebishev ? (((start + count) - (count - start) * std::cos( M_PI * (2 * i - 1) / (2 * count) )) / 2.)
                      : i
        );
    }

    return witnesses;
}

W_t::W_t(const witnesses_t &witnesses, const Circut &circut)
{
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
            m_map[witness] = condition;
        }
    }
}

W_t::cond_t W_t::operator()(witness_t w) const
{
    return *m_map.at(w);
}

GlobalParams::GlobalParams(const Circut &circut)
    //todo: генерация
    : m_TG{10, 20}
{
    m_witnesses = genWitnesses(wStart, circut.degree());

    generateT(circut);
    std::thread tS(&GlobalParams::generateS, this, std::ref(circut));
    std::thread tW(&GlobalParams::generateW, this, std::ref(circut));

    tS.join();
    tW.join();
}

witnesses_t GlobalParams::witnesses() const
{
    return m_witnesses;
}

witnesses_t GlobalParams::SWitnesses() const
{
    return m_SWitnesses;
}

TG_t GlobalParams::TG()
{
    return m_TG;
}

GlobalParams::ProverParams_t GlobalParams::PP()
{
    return {.TParams = {m_T, m_splittedT}, .SParams = {m_S, m_opsFromS}, .w = m_W};
}

void GlobalParams::generateT(const Circut &circut)
{
    dots_t dots, leftDots, rightDots, resultDots;

    auto cw = m_witnesses.cbegin();
    auto fillMap = [&cw, &dots](const values_t &row)
    {
        for(const auto& elment : row) {
            dots.push_back({X_t(*cw++), Y_t(elment)});
        }
    };

    fillMap(circut.m_inputX);
    fillMap(circut.m_inputW);

    m_SWitnesses = genWitnesses(wStart, circut.m_gates.size());
    auto i = m_SWitnesses.begin();

    for(const auto& gate : circut.m_gates) {
        dots.push_back({X_t(*cw++), Y_t(gate.m_input.a)});
        leftDots.push_back({X_t(*i), Y_t(gate.m_input.a)});

        dots.push_back({X_t(*cw++), Y_t(gate.m_input.b)});
        rightDots.push_back({X_t(*i), Y_t(gate.m_input.b)});

        dots.push_back({X_t(*cw++), Y_t(gate.m_output)});
        resultDots.push_back({X_t(*i), Y_t(gate.m_output)});

        i++;
    }

    m_T = T_t(dots);
    m_splittedT = {
                    .left = InterpolationPolynom(leftDots),
                    .right = InterpolationPolynom(rightDots),
                    .result = InterpolationPolynom(resultDots)
                  };
}

void GlobalParams::generateS(const Circut &circut)
{
    dots_t dots;

    for(std::size_t i = 0; i < circut.size(); i++) {
        auto currentOperation = circut.m_gates[i].m_type;
        auto sw = m_SWitnesses[i];

        dots.push_back({sw, currentOperation});

        FOROPS {
            m_opsFromS[operation].push_back(
                        operation == currentOperation ? dot_t{sw, 1} : dot_t{sw, 0}
            );
        }
    }

    m_S = S_t(dots);
}

void GlobalParams::generateW(const Circut &circut)
{
    m_W = W_t(m_witnesses, circut);
}

}
