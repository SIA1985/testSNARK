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

GlobalParams::GlobalParams(const Circut &circut)
    //todo: генерация
    : m_TG{10, 20}
{
    m_witnesses = genWitnesses(wStart, circut.degree(), wStep);

    auto mappedT = generateT(circut);
    std::thread tS(&GlobalParams::generateS, this, std::ref(circut));
    std::thread tW(&GlobalParams::generateW, this, std::ref(circut), std::ref(mappedT));

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
    return {.TParams = {m_T, m_splittedT}, .SParams = {m_S, m_opsFromS}, .WParams = {m_W, m_WT}};
}

GlobalParams::WtoValue_t GlobalParams::generateT(const Circut &circut)
{
    dots_t dots, leftDots, rightDots, resultDots;
    WtoValue_t mappedT;

    auto cw = m_witnesses.cbegin();
    auto fillMap = [&cw, &dots, &mappedT](const values_t &row)
    {
        for(const auto& element : row) {
            dots.push_back({*cw, element});
            mappedT[*cw] = element;
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
        mappedT[*cw] = gate.m_input.a;
        cw++;

        dots.push_back({*cw, gate.m_input.b});
        rightDots.push_back({*i, gate.m_input.b});
        mappedT[*cw] = gate.m_input.b;
        cw++;

        dots.push_back({*cw, gate.m_output});
        resultDots.push_back({*i, gate.m_output});
        mappedT[*cw] = gate.m_output;
        cw++;

        i++;
    }

    m_T = T_t(dots);
    m_splittedT = {
                    .left = InterpolationPolynom(leftDots),
                    .right = InterpolationPolynom(rightDots),
                    .result = InterpolationPolynom(resultDots)
                  };

    return mappedT;
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

void GlobalParams::generateW(const Circut &circut, WtoValue_t mappedT)
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

    dots_t dots;
    dots_t dotsWT;
    dots.reserve(duplicates.size());
    dotsWT.reserve(duplicates.size());

    for(const auto &[_, condition] : duplicates) {
        auto begin = condition->begin();
        auto end = condition->end();

        for(auto it = begin; it != end;) {
            dotsWT.push_back({*it, mappedT[*it]});
            dots.push_back({*it, *circleIterator(begin, ++it, end)});
        }
    }

    m_W = W_t(dots);
    m_WT = WT_t(dotsWT);

}

}
