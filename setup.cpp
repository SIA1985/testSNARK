#include "setup.h"
#include "operation.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <thread>

namespace snrk {

witnesses_t genWitnesses(witness_t start, std::size_t count, witness_t wStep)
{
    witnesses_t witnesses;
    for(witness_t i = 0; i < count; i++) {
        witnesses.push_back(start);
        start += wStep;
    }

    return witnesses;
}

GlobalParams::GlobalParams(const Circut &circut, const GPK_t &GPK)
    : m_GPK{GPK}
    , m_circutSize{circut.size()}
{
    m_witnesses = genWitnesses(wStart, circut.degree(), wStep);

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

std::size_t GlobalParams::witnessesCount() const
{
    return m_witnesses.size();
}

GPK_t GlobalParams::GPK() const
{
    return m_GPK;
}

GlobalParams::ProverParams_t GlobalParams::PP() const
{
    return {.TParams = {m_T, m_splittedT}, .SParams = {m_S, m_opsFromS}, .WParams = {m_WT, m_WI}};
}

std::size_t GlobalParams::circutSize() const
{
    return m_circutSize;
}

void GlobalParams::generateT(const Circut &circut)
{
    dots_t dots, leftDots, rightDots, resultDots;

    auto cw = m_witnesses.cbegin();
    auto fillMap = [&cw, &dots](const values_t &row)
    {
        for(const auto& element : row) {
            dots.push_back({*cw, element});
            cw++;
        }
    };

    fillMap(circut.m_inputX);
    fillMap(circut.m_inputW);

    m_SWitnesses = genWitnesses(wStart, circut.m_gates.size(), wStep);
    auto i = m_SWitnesses.begin();

    for(const auto& gate : circut.m_gates) {
        dots.push_back({*cw, gate.m_input.a});
        leftDots.push_back({*i, gate.m_input.a});
        cw++;

        dots.push_back({*cw, gate.m_input.b});
        rightDots.push_back({*i, gate.m_input.b});
        cw++;

        dots.push_back({*cw, gate.m_output});
        resultDots.push_back({*i, gate.m_output});
        cw++;

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

        for(auto operation : {Sum, Product}) {
            m_opsFromS[operation].push_back(
                        operation == currentOperation ? dot_t{sw, 1} : dot_t{sw, 0}
            );
        }
    }

    m_S = S_t(dots);
}

void GlobalParams::generateW(const Circut &circut)
{
    /*addr -> свидетели*/
    std::unordered_map<std::size_t, std::shared_ptr<cond_t>> duplicates;
    auto insert = [&duplicates](value_t key, witness_t value)
    {
        auto address = (std::size_t)key.address();
        if (!duplicates[address]) {
            duplicates[address] = std::make_shared<cond_t>();
        }

        duplicates[address]->insert(value);
    };

    auto cw = m_witnesses.cbegin();

    auto fillMap = [&cw, &insert](const values_t &row)
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


    auto circleIterator = [](cond_t::iterator begin, cond_t::iterator it, cond_t::iterator end)
    {
        if (it == end) {
            return begin;
        }

        return it;
    };

    dots_t dotsWI;
    dots_t dotsWT;
    dotsWI.reserve(duplicates.size());
    dotsWT.reserve(duplicates.size());

    for(const auto &[_, condition] : duplicates) {
        auto begin = condition->begin();
        auto end = condition->end();

        for(auto it = begin; it != end;) {
            dotsWT.push_back({*it, *it});
            dotsWI.push_back({*it, *circleIterator(begin, ++it, end)});
        }
    }

    m_WT = W_t(dotsWT); //witness -> witness
    m_WI = W_t(dotsWI); //witness -> next_this_value_wintess
}

}
